import copy
import os

def build_params(**module_templates):

     env = os.environ
     schedule_file=f"schedules/a2a_qed_and_vec_meson_{env['SERIES']}{env['CFG']}.sched"
     
     sources = str(3*int(env['NOISE'])*int(env["TIME"]))

     gammas = {
          #"pion_local"  :["(G5 G5)"],
          "vec_local"   :" ".join(["(GX GX)","(GY GY)","(GZ GZ)"]),
          #"vec_onelink" :["(GX G1)","(GY G1)","(GZ G1)"]
     }
     gamma_label, gamma_string = next(iter(gammas.items()))
     
     currents = {
         "current_local"  :" ".join(["(GX GX)","(GY GY)","(GZ GZ)","(GT GT)"]),
         #"current_onelink": " ".join(["(GX G1)","(GY G1)","(GZ G1)","(GT G1)"])
     }
     current_label, current_string = next(iter(currents.items()))

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
                 "schedule_file":schedule_file,
                 "saveSchedule":"false",
                 "parallelWriteMaxRetry":"-1",
             },
             "modules":{},
         },
     }

     modules = []

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge_fat"
     module["options"]["file"] = f'configs/fat{env["ENS"]}{env["SERIES"]}.ildg'
     modules.append(module)

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge_long"
     module["options"]["file"] = f'configs/lng{env["ENS"]}{env["SERIES"]}.ildg'
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

     module = copy.deepcopy(module_templates["lma_solver"])
     module["id"]["name"] = f"project_low"
     module["options"]["action"] = f"stag_{mass_string}"
     module["options"]["projector"] = "true"
     module["options"]["lowModes"] = f"evecs_{mass_string}"
     modules.append(module)

     module = copy.deepcopy(module_templates["em_func"])
     module["id"]["name"] = f"photon_func"
     module["options"]["gauge"] = "feynman"
     module["options"]["zmScheme"] = "qedL"
     modules.append(module)

     seed=env["SEEDSTRING"]
     for i,sw in enumerate(range(int(env["WSEEDSTART"]),int(env["WSEEDSTART"])+int(env["NSEEDS"]))):

          wseed=f"w{seed}{sw}"
          module = copy.deepcopy(module_templates["full_volume_noise"])
          module["id"]["name"] = f"noise_{seed}{sw}"
          module["options"]["nsrc"] = env['NOISE']
          modules.append(module)

          for j,sv in enumerate(range(int(env["VSEEDSTART"]),int(env["VSEEDSTART"])+int(env["NSEEDS"]))):

               vseed=f"v{seed}{sv}"
               if i == 0:
                    module = copy.deepcopy(module_templates["load_vectors"])
                    module["id"]["name"] = vseed
                    module["options"]["filestem"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/vectors/{mass_string}/{seed}{sv}_v"
                    module["options"]["multiFile"] = 'false'
                    module['options']['size'] = sources
                    modules.append(module)
          
               if sv == sw:
                    continue
               
               module = copy.deepcopy(module_templates["qed_meson_field"])
               module["id"]["name"] = f"mf_{wseed}_{vseed}_qed"
               module["options"].update({
                    "action":f"stag_{mass_string}",
                    "block":"500",
                    "left":f"noise_{seed}{sw}_vec",
                    "right":vseed,
                    "em_func":"photon_func",
                    "EmSeedString":env['EMSEEDSTRING'],
                    "nem_fields":env['NEM'],
                    "spinTaste":{
                         "gammas":current_string,
                         "gauge" :"",
                         "applyG5": "false"
                    },
                    "lowModes":f"",
                    "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/test/mesons/a2a/mf_{env['SERIES']}_{wseed}_{vseed}"
               })
               modules.append(module)                    
         
               module = copy.deepcopy(module_templates["meson_field"])
               module["id"]["name"] = f"mf_{wseed}_{vseed}"
               module["options"].update({
                    "action":f"stag_{mass_string}",
                    "block":"500",
                    "left":f"noise_{seed}{sw}_vec",
                    "right":vseed,
                    "spinTaste":{
                         "gammas":gamma_string,
                         "gauge" :"gauge" if 'onelink' in gamma_label else  "",
                         "applyG5":"false"
                    },
                    "lowModes":f"",
                    "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/test/mesons/a2a/mf_{env['SERIES']}_{wseed}_{vseed}"
               })
               modules.append(module)                    

     params["grid"]["modules"] = {"module":modules}
     
     moduleList = [m["id"]["name"] for m in modules]

     f = open(schedule_file, "w")
     f.write(str(len(moduleList)) + "\n" + "\n".join(moduleList))
     f.close()

     return params
