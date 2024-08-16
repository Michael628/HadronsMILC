import copy
import os
import random

def buildParams(**moduleTemplates):

    env = os.environ

    jobid=int(random.random()*100)
    scheduleFile=f"schedules/lma_{jobid}.sched"
    masses=env["MASSES"].strip().split(" ")
    gammas = {
        "pion_local"   :"(G5 G5)",
        "vec_local"  :" ".join(["(GX GX)","(GY GY)","(GZ GZ)"])
    }
    gammas_iter = list(gammas.items())

    # Make sure we iterate over pion first
    gammas_iter.sort(key=(lambda a: a[0] != "pion_local"))
    
    params = {
        "grid":{
            "parameters":{
                "runId":f"LMI-RW-series-{env['SERIES']}-{env['EIGS']}-eigs-{env['NOISE']}-noise",
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
                f"scheduleFile":scheduleFile,
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
             "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f"configs/fat{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f"configs/lng{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["castGauge"])
    module["id"]["name"] = "gauge_fatf"
    module["options"]["field"] = "gauge_fat"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["castGauge"])
    module["id"]["name"] = "gauge_longf"
    module["options"]["field"] = "gauge_long"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["action"])
    module["id"]["name"] = "stag_e"
    module["options"]["mass"] = "0.0"
    module["options"]["gaugefat"] = "gauge_fat"
    module["options"]["gaugelong"] = "gauge_long"
    modules.append(module)

     
    module = copy.deepcopy(moduleTemplates["op"])
    module["id"]["name"] = "stag_op"
    module["options"]["action"] = "stag_e"
    modules.append(module)

    alpha = str(float( "%0.2g" % (float(env['ALPHA'])+float(env["BUNDLEINDEX"])*float(env["DALPHA"])))) if "DALPHA" in env.keys() else env['ALPHA']
    npoly = str(int(env['NPOLY'])+int(env["BUNDLEINDEX"])*int(env["DNPOLY"])) if "DNPOLY" in env.keys() else env['NPOLY']

    module = copy.deepcopy(moduleTemplates["IRL"])
    module["id"]["name"] = "epack"
    module["id"]["type"] = "MSolver::StagIRL"
    module["options"]["op"] = "stag_op_schur"
    module["options"]["lanczosParams"]["Cheby"]["alpha"] = alpha
    module["options"]["lanczosParams"]["Cheby"]["beta"] = env['BETA']
    module["options"]["lanczosParams"]["Cheby"]["Npoly"] = npoly
    module["options"]["lanczosParams"]["Nstop"] = env['NSTOP']
    module["options"]["lanczosParams"]["Nk"] = env['NK']
    module["options"]["lanczosParams"]["Nm"] = env['NM']
    if 'EIGOUT' in env.keys():
        module["options"]["output"] = env['EIGOUT'].format(ens=env['ENS'],series=env['SERIES'])
    if 'MULTIFILE' in env.keys():
        module["options"]["multiFile"] = env['MULTIFILE']
    if 'EIGRESID' in env.keys():
        module["options"]["lanczosParams"]["resid"] = env['EIGRESID']
    modules.append(module)

    time = int(env["TIME"])
    tStart = int(env["TSTART"])
    tStop = int(env["TSTOP"])+1
    tStep = int(env["DT"])
    for time_index in range(tStart,tStop,tStep):
        block_label=f"t{time_index}"
        noise=f"noise_{block_label}"

        module = copy.deepcopy(moduleTemplates["noiseRW"])
        module["id"]["name"] = noise
        module["options"]["nSrc"] = env["NOISE"]
        module["options"]["tStep"] = str(time)
        module["options"]["t0"] = str(time_index)
        modules.append(module)      

        for mass1_index, m1 in enumerate(masses):
            mass1 = "0." + m1
            mass1_label = "m"+m1

            if time_index == tStart:
                module = copy.deepcopy(moduleTemplates["epackModify"])
                module["id"]["name"] = f"evecs_{mass1_label}"
                module["options"]["eigenPack"] = "epack"
                module["options"]["mass"] = mass1
                modules.append(module)
                 
                module = copy.deepcopy(moduleTemplates["action"])
                module["id"]["name"] = f"stag_{mass1_label}"
                module["options"]["mass"] = mass1
                module["options"]["gaugefat"] = "gauge_fat"
                module["options"]["gaugelong"] = "gauge_long"
                modules.append(module)

                if mass1_index == 0:
                    module = copy.deepcopy(moduleTemplates["mesonField"])
                    module["id"]["name"] = f"mf_ll_wv_local"
                    module["options"].update({
                        "action":f"stag_{mass1_label}",
                        "block":"200",
                        "spinTaste":{
                            "gammas":"(G5 G5) (GX GX) (GY GY) (GZ GZ)",
                            "gauge" :"",
                            "applyG5":"false"
                        },
                        "lowModes":f"evecs_{mass1_label}",
                        "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{(mass1_label+'/') if len(masses) > 1 else ''}mf_{env['SERIES']}"
                    })
                    modules.append(module)

                module = copy.deepcopy(moduleTemplates["actionF"])
                module["id"]["name"] = f"istag_{mass1_label}"
                module["options"]["mass"] = mass1
                module["options"]["gaugefat"] = "gauge_fatf"
                module["options"]["gaugelong"] = "gauge_longf"
                modules.append(module)
         
                module = copy.deepcopy(moduleTemplates["cgMP"])
                module["id"]["name"] = f"stag_ama_{mass1_label}"
                module["options"]["outerAction"] = f"stag_{mass1_label}"
                module["options"]["innerAction"] = f"istag_{mass1_label}"
                module["options"]["residual"] = "1e-8"
                modules.append(module)
          
                module = copy.deepcopy(moduleTemplates["lmaProj"])
                module["id"]["name"] = f"stag_ranLL_8000_{mass1_label}"
                module["options"]["action"] = f"stag_{mass1_label}"
                module["options"]["lowModes"] = f"evecs_{mass1_label}"
                module["options"]["eigStart"] = "0"
                module["options"]["nEigs"] = "8000"
                modules.append(module)

                module = copy.deepcopy(moduleTemplates["lmaProj"])
                module["id"]["name"] = f"stag_ranLL_4000_{mass1_label}"
                module["options"]["action"] = f"stag_{mass1_label}"
                module["options"]["lowModes"] = f"evecs_{mass1_label}"
                module["options"]["eigStart"] = "0"
                module["options"]["nEigs"] = "4000"
                modules.append(module)

                module = copy.deepcopy(moduleTemplates["lmaProj"])
                module["id"]["name"] = f"stag_ranLL_2000_{mass1_label}"
                module["options"]["action"] = f"stag_{mass1_label}"
                module["options"]["lowModes"] = f"evecs_{mass1_label}"
                module["options"]["eigStart"] = "0"
                module["options"]["nEigs"] = "2000"
                modules.append(module)

            for gamma_label, gamma_string in gammas_iter:

                # Only do sea mass for onelink
                if "onelink" in gamma_label and mass1_index != 0:
                    continue

                # Attach gauge link field to spin-taste for onelink
                gauge = "gauge" if "onelink" in gamma_label else ""

                for solver_label in ["ranLL_2000","ranLL_4000","ranLL_8000","ama"]:

                    solver=f"stag_{solver_label}_{mass1_label}"
                    guess=f"quark_ranLL_8000_{gamma_label}_{mass1_label}_{block_label}" if solver_label == "ama" else ""
                    quark_m1=f"quark_{solver_label}_{gamma_label}_{mass1_label}_{block_label}"
                    
                    module = copy.deepcopy(moduleTemplates["quarkProp"])
                    module["id"]["name"] = quark_m1
                    module["options"].update({
                        "source"   :noise,
                        "solver"   :solver,
                        "guess"    :guess,
                        "spinTaste":{
                            "gammas":gamma_string,
                            "gauge" :gauge,
                            "applyG5":"true"
                        }
                    })
                    modules.append(module)

                    # Perform all cross-contractions between masses
                    for mass2_index, m2 in enumerate(masses):

                        if mass2_index > mass1_index:
                            continue

                        mass2 = "0." + m2
                        mass2_label="m"+m2

                        if mass2_index == mass1_index:
                            mass_label = mass1_label
                        else:
                            mass_label = mass1_label + "_" + mass2_label

                        quark_m2=f"quark_{solver_label}_pion_local_{mass2_label}_{block_label}"

                        module = copy.deepcopy(moduleTemplates["propContract"])
                        module["id"]["name"] = f"corr_{solver_label}_{gamma_label}_{mass_label}_{block_label}"
                        module["options"].update({
                            "source":quark_m1,
                            "sink":quark_m2,
                            "sinkFunc":"sink",
                            "sourceShift":noise+"_shift",
                            "sourceGammas":gamma_string,
                            "sinkSpinTaste":{
                                "gammas":gamma_string,
                                "gauge" :gauge,
                                "applyG5":"true"
                            },
                            "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/correlators/{mass_label+'/' if len(masses) > 1 else ''}{gamma_label}/{solver_label}/corr_{gamma_label}_{solver_label}_{mass_label}_{block_label}_{env['SERIES']}",
                        })
                        modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]
    moduleList.sort(key=lambda x: 1 if "mf" in x else 0)
    
    f = open(scheduleFile,"w")
    f.write(str(len(moduleList))+"\n"+"\n".join(moduleList))
    f.close()
    
    return params
