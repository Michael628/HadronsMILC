import copy
import os

def buildParams(**moduleTemplates):

    env = os.environ
    scheduleFile=f"schedules/eig_lma1_meson_{env['SERIES']}{env['CFG']}.sched"

    masses=env["MASSES"].strip().split(" ")

    gammas = {
         "pion_local"   :"(G5 G5)",
         "vec_local"  :" ".join(["(GX GX)","(GY GY)","(GZ GZ)"]),
         "vec_onelink": " ".join(["(GX G1)","(GY G1)","(GZ G1)"])
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
                "scheduleFile":scheduleFile,
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
            "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge"
    module["options"]["file"] = f"lat/scidac/l{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)
    
    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f"lat/scidac/fat{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f"lat/scidac/lng{env['ENS']}{env['SERIES']}.ildg"
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

    module = copy.deepcopy(moduleTemplates["epackLoad"])
    module["id"]["name"] = "epack2k"
    module["options"]["filestem"] = f"eigen/eig{env['ENS']}nv2000er8_grid_{env['SERIES']}"
    module["options"]["size"] = env['SOURCEEIGS']
    module["options"]["multiFile"] = "true"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["epackModify"])
    module["id"]["name"] = f"evecs_m0"
    module["options"]["eigenPack"] = "epack2k"
    module["options"]["mass"] = "0.0"
    modules.append(module)
                 
    module = copy.deepcopy(moduleTemplates["lmaProj"])
    module["id"]["name"] = f"stag_ranLL_epack"
    module["options"]["action"] = f"stag_e"
    module["options"]["lowModes"] = f"evecs_m0"
    module["options"]["projector"] = "true"
    module["options"]["eigStart"] = "0"
    module["options"]["nEigs"] = env['SOURCEEIGS']
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["IRL"])
    module["id"]["name"] = "epack"
    module["options"]["op"] = "stag_op_schur"
    module["options"]["lanczosParams"]["Cheby"]["alpha"] = env['ALPHA']
    module["options"]["lanczosParams"]["Cheby"]["beta"] = env['BETA']
    module["options"]["lanczosParams"]["Cheby"]["Npoly"] = env['NPOLY']
    module["options"]["lanczosParams"]["Nstop"] = env['NSTOP']
    module["options"]["lanczosParams"]["Nk"] = env['NK']
    module["options"]["lanczosParams"]["Nm"] = env['NM']
    module["options"]["lanczosParams"]["resid"] = env['RESID'] if 'RESID' in env else '1e-8'
    module["options"]["epackIn"] = 'epack2k'
    module["options"]["projector"] = 'stag_ranLL_epack'
    module["options"]["output"] = env['EIGOUT'].format(ens=env['ENS'],series=env['SERIES'])
    module["options"]["multiFile"] = "true"
    modules.append(module)
    
    module = copy.deepcopy(moduleTemplates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
    modules.append(module)

    time = int(env["TIME"])
    tStart = int(env["TSTART1"])
    tStop = int(env["TSTOP1"])+1
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
                    module["id"]["name"] = f"mf_ll_wv_onelink"
                    module["options"].update({
                        "action":f"stag_{mass1_label}",
                        "block":env['BLOCKSIZE'],
                        "spinTaste":{
                            "gammas":"(GX G1) (GY G1) (GZ G1)",
                            "gauge" :"gauge",
                            "applyG5":"false"
                        },
                        "lowModes":f"evecs_{mass1_label}",
                        "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{(mass1_label+'/') if len(masses) > 1 else ''}mf_{env['SERIES']}"
                    })
                    modules.append(module)

                    module = copy.deepcopy(moduleTemplates["mesonField"])
                    module["id"]["name"] = f"mf_ll_wv_local"
                    module["options"].update({
                        "action":f"stag_{mass1_label}",
                        "block":env['BLOCKSIZE'],
                        "spinTaste":{
                            "gammas":"(G5 G5) (GX GX) (GY GY) (GZ GZ)",
                            "gauge" :"",
                            "applyG5":"false"
                        },
                        "lowModes":f"evecs_{mass1_label}",
                        "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{(mass1_label+'/') if len(masses) > 1 else ''}mf_{env['SERIES']}"
                    })
                    modules.append(module)

                module = copy.deepcopy(moduleTemplates["lmaProj"])
                module["id"]["name"] = f"stag_ranLL_{mass1_label}"
                module["options"]["action"] = f"stag_{mass1_label}"
                module["options"]["lowModes"] = f"evecs_{mass1_label}"
                module["options"]["eigStart"] = "0"
                module["options"]["nEigs"] = env['EIGS']
                modules.append(module)

            for gamma_label, gamma_string in gammas_iter:

                # Only do sea mass for onelink
                if "onelink" in gamma_label and mass1_index != 0:
                    continue

                # Attach gauge link field to spin-taste for onelink
                gauge = "gauge" if "onelink" in gamma_label else ""

                for solver_label in ["ranLL"]:

                    solver=f"stag_{solver_label}_{mass1_label}"
                    guess=f"quark_ranLL_{gamma_label}_{mass1_label}_{block_label}" if solver_label == "ama" else ""
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
                            "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/correlators/{mass_label}/{gamma_label}/{solver_label}/corr_{gamma_label}_{solver_label}_{mass_label}_{block_label}_{env['SERIES']}",
                        })
                        modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]

    f = open(scheduleFile,"w")
    f.write(str(len(moduleList))+"\n"+"\n".join(moduleList))
    f.close()
     
    return params
