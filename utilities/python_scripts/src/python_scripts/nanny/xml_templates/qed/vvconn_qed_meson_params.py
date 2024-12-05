import copy
import os

def build_params(**module_templates):

     env = os.environ
     schedule_file=f"schedules/vvconn_qed_meson_{env['SERIES']}{env['CFG']}.sched"
     
     gammas = {
          "pion_local"  :["(G5 G5)"],
          "vec_local"   :["(GX GX)","(GY GY)","(GZ GZ)"],
          #"vec_onelink" :["(GX G1)","(GY G1)","(GZ G1)"]
     }
     # The key to use when contracting <pion solve| vec_gamma | a2a solve>
     vec_gamma_key = list(filter(lambda x: True if "vec" in x else False,gammas.keys()))[0]
     
     gammas_iter = list(gammas.items())

     currents = {
         "current_local"  :" ".join(["(GX GX)","(GY GY)","(GZ GZ)","(GT GT)"]),
         #"current_onelink": " ".join(["(GX G1)","(GY G1)","(GZ G1)","(GT G1)"])
     }
     
     # Make sure we iterate over pion first
     gammas_iter.sort(key=(lambda a: a[0] != "pion_local"))

     sources = str(3*int(env['NOISE'])*int(env["TIME"]))

     params = {
         "grid":{
             "parameters":{
                 "runId":f"A2A-series-{env['SERIES']}-{env['EIGS']}-eigs-{env['NOISE']}-noise",
                 "trajCounter":{
                     "start":env["CFG"],
                     "end":"10000",
                     "step":"10000",
                 },
                 "genetic":{
                     "popSize":"20",
                     "maxGen":"1000",
                     "maxCstGen":"100",
                     "mutationRate":"0.1",
                 },
                 "graphFile":"",
                 "scheduleFile":schedule_file,
                 "saveSchedule":"false",
                 "parallelWriteMaxRetry":"-1",
             },
             "modules":{},
         },
     }

     modules = []

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge_fat"
     module["options"]["file"] = f'lat/scidac/fat{env["ENS"]}{env["SERIES"]}.ildg'
     modules.append(module)

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge_long"
     module["options"]["file"] = f'lat/scidac/lng{env["ENS"]}{env["SERIES"]}.ildg'
     modules.append(module)

     module = copy.deepcopy(module_templates["cast_gauge"])
     module["id"]["name"] = "gauge_fatf"
     module["options"]["field"] = "gauge_fat"
     modules.append(module)

     module = copy.deepcopy(module_templates["cast_gauge"])
     module["id"]["name"] = "gauge_longf"
     module["options"]["field"] = "gauge_long"
     modules.append(module)

     module = copy.deepcopy(module_templates["epack_load"])
     module["id"]["name"] = "epack"
     module["options"]["filestem"] = f'eigs/eig{env["ENS"]}nv{env["SOURCEEIGS"]}{env["SERIES"]}'
     module["options"]["size"] = env["EIGS"]
     modules.append(module)

     m = os.environ["MASSES"].strip().split(" ")[0]
     mass_string = f"m{m}"
     mass = f"0.{m}"

     module = copy.deepcopy(module_templates["action"])
     module["id"]["name"] = f"stag_{mass_string}"
     module["options"]["mass"] = mass
     module["options"]["gaugefat"] = "gauge_fat"
     module["options"]["gaugelong"] = "gauge_long"
     modules.append(module)

     module = copy.deepcopy(module_templates["epack_modify"])
     module["id"]["name"] = f"evecs_{mass_string}"
     module["options"]["eigenPack"] = "epack"
     module["options"]["mass"] = mass
     modules.append(module)

     module = copy.deepcopy(module_templates["action_float"])
     module["id"]["name"] = f"istag_{mass_string}"
     module["options"]["mass"] = mass
     module["options"]["gaugefat"] = "gauge_fatf"
     module["options"]["gaugelong"] = "gauge_longf"
     modules.append(module)
         
     module = copy.deepcopy(module_templates["lma_solver"])
     module["id"]["name"] = f"stag_ranLL_{mass_string}"
     module["options"]["action"] = f"stag_{mass_string}"
     module["options"]["lowModes"] = f"evecs_{mass_string}"
     modules.append(module)

     module = copy.deepcopy(module_templates["mixed_precision_cg"])
     module["id"]["name"] = f"stag_ama_{mass_string}"
     module["options"]["outerAction"] = f"stag_{mass_string}"
     module["options"]["innerAction"] = f"istag_{mass_string}"
     module["options"]["residual"] = "1e-8"
     modules.append(module)
          
     module = copy.deepcopy(module_templates["em_func"])
     module["id"]["name"] = f"photon_func"
     module["options"]["gauge"] = "feynman"
     module["options"]["zmScheme"] = "qedL"
     modules.append(module)

     seed=env["SEEDSTRING"]
     ext_seed = f"ext_{seed}0"

     gamma_label = 'pion_local'
     for gamma in gammas[gamma_label]:
          gamma_suf = gamma.replace("(","")
          gamma_suf = gamma_suf.replace(")","")
          gamma_suf = gamma_suf.replace(" ","_")

          module = copy.deepcopy(module_templates["load_vectors"])
          module["id"]["name"] = f"{ext_seed}_{gamma_suf}"
          module["options"]["filestem"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/ext_vectors/{mass_string}/{gamma_suf}/{ext_seed}_v"
          module["options"]["multiFile"] = 'false'
          module['options']['size'] = str(3*int(env['NOISE'])*int(env["TIME"]))
          modules.append(module)
          
          for i,s in enumerate(range(int(env["NSEEDS"]))):

               vseed = f"v{seed}{s}"
               module = copy.deepcopy(module_templates["load_vectors"])
               module["id"]["name"] = vseed
               module["options"]["filestem"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/vectors/{mass_string}/{seed}{s}_v"
               module["options"]["multiFile"] = 'false'
               module['options']['size'] = str(3*int(env['NOISE'])*int(env["TIME"]))
               modules.append(module)
          
               # Builds meson fields of dimension (Nt,NSEEDS,Nt). To time dilute, the first two dims become the diagonal
               # of a meson field with dim (Nt,Nt*NSEEDS,Nt).
               module = copy.deepcopy(module_templates["qed_meson_field"])
               module["id"]["name"] = f"{gamma_suf}_{ext_seed}_G5herm_em_{vseed}"
               module["options"].update({
                    "action":f"stag_{mass_string}",
                    "block":"500",
                    "left":f"{ext_seed}_{gamma_suf}",
                    "right":vseed,
                    "em_func":"photon_func",
                    "EmSeedString":env['EMSEEDSTRING'],
                    "nem_fields":env['NEM'],
                    "spinTaste":{
                         "gammas":currents['current_local'],
                         "gauge" :"",
                         "applyG5": "true"
                    },
                    "lowModes":"",
                    "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/test/mesons/ext_G5herm/mf_{gamma_suf}{ext_seed}_G5herm_{vseed}_{env['SERIES']}"
               })
               modules.append(module)
                     
               module = copy.deepcopy(module_templates["meson_field"])
               module["id"]["name"] = f"{gamma_suf}_{ext_seed}_G5herm_vec_{vseed}"
               module["options"].update({
                    "action":f"stag_{mass_string}",
                    "block":"500",
                    "left":f"{ext_seed}_{gamma_suf}",
                    "right":vseed,
                    "spinTaste":{
                         "gammas":" ".join(gammas[vec_gamma_key]),
                         "gauge" :"gauge" if "onelink" in vec_gamma_key else "",
                         "applyG5":"true"
                    },
                    "lowModes":f"",
                    "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/test/mesons/ext_G5herm/mf_{gamma_suf}{ext_seed}_G5herm_{vseed}_{env['SERIES']}"
               })
               modules.append(module)
               
     gamma_label = vec_gamma_key
     for i,gamma in enumerate(gammas[gamma_label]):
          gamma_suf = gamma.replace("(","")
          gamma_suf = gamma_suf.replace(")","")
          gamma_suf = gamma_suf.replace(" ","_")

          module = copy.deepcopy(module_templates["load_vectors"])
          module["id"]["name"] = f"{ext_seed}_{gamma_suf}"
          module["options"]["filestem"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/ext_vectors/{mass_string}/{gamma_suf}/{ext_seed}_v"
          module["options"]["multiFile"] = 'false'
          module['options']['size'] = str(3*int(env['NOISE'])*int(env["TIME"]))
          modules.append(module)
          
          for s in range(int(env["NSEEDS"])):

               wseed = f"w{seed}{s}"
               if i == 0:
                    # Build noise sources without time dilution
                    module = copy.deepcopy(module_templates["full_volume_noise"])
                    module["id"]["name"] = f"noise_{seed}{s}"
                    module["options"]["nsrc"] = env["NOISE"]
                    modules.append(module)

               # Builds meson fields of dimension (Nt,NSEEDS,Nt). To time dilute, the first two dims become the diagonal
               # of a meson field with dim (Nt,Nt*NSEEDS,Nt).
               module = copy.deepcopy(module_templates["qed_meson_field"])
               module["id"]["name"] = f"mf_{wseed}_{ext_seed}_{gamma_suf}"
               module["options"].update({
                    "action":f"stag_{mass_string}",
                    "block":"500",
                    "left":f"noise_{seed}{s}_vec",
                    "right":f"{ext_seed}_{gamma_suf}",
                    "em_func":"photon_func",
                    "EmSeedString":env['EMSEEDSTRING'],
                    "nem_fields":env['NEM'],
                    "spinTaste":{
                         "gammas":currents['current_local'],
                         "gauge" :"",
                         "applyG5": "false"
                    },
                    "lowModes":"",
                    "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/test/mesons/ext/mf_{wseed}_{ext_seed}{gamma_suf}_{env['SERIES']}"
               })
               modules.append(module)
         
     params["grid"]["modules"] = {"module":modules}
     
     moduleList = [m["id"]["name"] for m in modules]

     f = open(schedule_file, "w")
     f.write(str(len(moduleList)) + "\n" + "\n".join(moduleList))
     f.close()

     return params
