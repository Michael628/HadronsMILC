import copy
import os

def build_params(**module_templates):

     env = os.environ
     schedule_file=f"schedules/a2a_proj_meson_{env['SERIES']}{env['CFG']}.sched"
     
     sources = str(3*int(env['NOISE'])*int(env["TIME"]))

     gammas = ["(G1 G1)", "(G5 G5)", "(GX GX)","(GY GY)","(GZ GZ)"]
     gamma_string=" ".join(gammas)

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

     seed=env["SEEDSTRING"]
     for i,sw in enumerate(range(int(env["WSEEDSTART"]),int(env["WSEEDSTART"])+int(env["NSEEDS"]))):

          wseed=f"w{seed}{sw}"
          if env['DOHIGH'] == 'true':
               module = copy.deepcopy(module_templates["time_diluted_noise"])
               module["id"]["name"] = f"noise_{seed}{sw}"
               module["options"]["nsrc"] = env['NOISE']
               modules.append(module)

               module = copy.deepcopy(module_templates["quark_prop"])
               module["id"]["name"] = wseed
               module["options"].update({
                    "source"   :f"noise_{seed}{sw}_vec",
                    "solver"   :f"project_low_subtract",
                    "guess"    :"",
                    "spinTaste":{
                         "gammas":"",
                         "gauge" :"",
                         "applyG5":"false"
                    }
               })
               modules.append(module)
          
               if env['DOEIG'] == 'true':
                    module = copy.deepcopy(module_templates["meson_field"])
                    module["id"]["name"] = f"mf_{wseed}_eig"
                    module["options"].update({
                         "action":f"stag_{mass_string}",
                         "block":"500",
                         "left":wseed,
                         "right":"",
                         "spinTaste":{
                              "gammas":gamma_string,
                              "gauge" :"gauge" if "G1" in gamma_string and "G1 G1" not in gamma_string else "",
                              "applyG5":"false"
                         },
                         "lowModes":f"evecs_{mass_string}",
                         "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/a2a/mf_{env['SERIES']}_{wseed}_e{env['EIGS']}"
                    })
                    modules.append(module)                    

          for j,sv in enumerate(range(int(env["VSEEDSTART"]),int(env["VSEEDSTART"])+int(env["NSEEDS"]))):

               vseed=f"v{seed}{sv}"
               if i == 0:
                    if env['DOHIGH'] == 'true':
                         module = copy.deepcopy(module_templates["load_vectors"])
                         module["id"]["name"] = vseed
                         module["options"]["filestem"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/vectors/{mass_string}/{seed}{sv}_v"
                         module["options"]["multiFile"] = 'true'
                         module['options']['size'] = sources
                         modules.append(module)
          

                         if env['DOEIG'] == 'true':
                              module = copy.deepcopy(module_templates["meson_field"])
                              module["id"]["name"] = f"mf_eig_{vseed}"
                              module["options"].update({
                                   "action":f"stag_{mass_string}",
                                   "block":"500",
                                   "left":"",
                                   "right":vseed,
                                   "spinTaste":{
                                        "gammas":gamma_string,
                                        "gauge" :"gauge" if "G1" in gamma_string and "G1 G1" not in gamma_string else "",
                                        "applyG5":"false"
                                   },
                                   "lowModes":f"evecs_{mass_string}",
                                   "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/a2a/mf_{env['SERIES']}_e{env['EIGS']}_{vseed}"
                              })
                              modules.append(module)                    

                    if j == 0 and env['DOEIG'] == 'true':
                         module = copy.deepcopy(module_templates["meson_field"])
                         module["id"]["name"] = f"mf_eig_eig"
                         module["options"].update({
                              "action":f"stag_{mass_string}",
                              "block":"500",
                              "left":"",
                              "right":"",
                              "spinTaste":{
                                   "gammas":gamma_string,
                                   "gauge" :"gauge" if "G1" in gamma_string and "G1 G1" not in gamma_string else "",
                                   "applyG5":"false"
                              },
                              "lowModes":f"evecs_{mass_string}",
                              "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/a2a/mf_{env['SERIES']}_e{env['EIGS']}_e{env['EIGS']}"
                         })
                         modules.append(module)

               if sv == sw:
                    continue

               if env['DOHIGH'] == 'true':
                    module = copy.deepcopy(module_templates["meson_field"])
                    module["id"]["name"] = f"mf_{wseed}_{vseed}"
                    module["options"].update({
                         "action":f"stag_{mass_string}",
                         "block":"500",
                         "left":wseed,
                         "right":vseed,
                         "spinTaste":{
                              "gammas":gamma_string,
                              "gauge" :"gauge" if "G1" in gamma_string and "G1 G1" not in gamma_string else "",
                              "applyG5":"false"
                         },
                         "lowModes":f"",
                         "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/a2a/mf_{env['SERIES']}_{wseed}_{vseed}"
                    })
                    modules.append(module)                    
         
     params["grid"]["modules"] = {"module":modules}
     
     moduleList = [m["id"]["name"] for m in modules]

     f = open(schedule_file, "w")
     f.write(str(len(moduleList)) + "\n" + "\n".join(moduleList))
     f.close()

     return params
