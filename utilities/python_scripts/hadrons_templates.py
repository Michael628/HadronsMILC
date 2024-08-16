import copy

def generateTemplates():
    moduleTemplates = {}

    # Loads gauge fields from ILDG file
    moduleTemplates["loadGauge"] = {
        "id":{
            "name":"",
            "type": "MIO::LoadIldg"
        },
        "options": {
            "file":""
        }
    }
    
    # Recasts double precision gauge field to single
    moduleTemplates["castGauge"] = {
        "id":{
            "name":"",
            "type": "MUtilities::GaugeSinglePrecisionCast"
        },
        "options": {
            "field":""
        }
    }

    # Creates Dirac Matrix
    moduleTemplates["action"] = {
        "id":{
            "name":"",
            "type":"MAction::ImprovedStaggeredMILC"
        },
        "options":{
            "mass":"",
            "c1":"1.0",
            "c2":"1.0",
            "tad":"1.0",
            "boundary":"1 1 1 1",
            "twist":"0 0 0",
            "Ls":"1",
            "gaugefat":"",
            "gaugelong":""
        }
    }

    # Creates single precision Dirac Matrix
    moduleTemplates["actionF"] = copy.deepcopy(moduleTemplates["action"])
    moduleTemplates["actionF"]["id"]["type"] = "MAction::ImprovedStaggeredMILCF"

    moduleTemplates["op"] = {
        "id":{
            "name":"",
            "type":"MFermion::StagOperators"
        },
        "options":{
            "action":""
        }
    }

    moduleTemplates["IRL"] = {
        "id":{
            "name":"",
            "type":"MSolver::StagFermionIRL"
        },
        "options":{
            "op": "",
            "lanczosParams": {
                "Cheby": {
                    "alpha":"",
                    "beta":"",
                    "Npoly":""
                },
                "Nstop":"",
                "Nk":"",
                "Nm":"",
                "resid":"1e-8",
                "MaxIt":"5000",
                "betastp":"0",
                "MinRes":"0"
            },
            "evenEigen":"false",
            "redBlack":"true",
            "output":"",
            "multiFile":"false"
        }
    }
        
    # Loads eigenvector pack
    moduleTemplates["epackLoad"] = {
        "id":{
            "name":"",
            "type":"MIO::StagLoadFermionEigenPack"
        },
        "options":{
            "redBlack":"true",
            "filestem":"",
            "multiFile":"false",
            "size":"",
            "Ls":"1"
        }
    }

    # Write eigenvalues into separate XML file
    moduleTemplates["evalSave"] = {
        "id":{
            "name":"",
            "type":"MUtilities::EigenPackExtractEvals"
        },
        "options":{
            "eigenPack":"",
            "output":""
        }
    }

    # Shift eigenvalues of eigenpack by a mass value
    moduleTemplates["epackModify"] = {
        "id":{
            "name":"",
            "type":"MUtilities::ModifyEigenPackMILC"
        },
        "options":{
            "eigenPack":"",
            "mass":"",
            "evenEigen":"false",
            "normalizeCheckerboard":"false"
        }
    }

    # Create operator to apply spin taste to fields
    moduleTemplates["spinTaste"] = {
        "id":{
            "name":"",
            "type":"MFermion::SpinTaste"
        }
    }

    # Creates function that spacially ties up field with chosen momentum phase
    moduleTemplates["sink"] = {
        "id":{
            "name":"",
            "type":"MSink::ScalarPoint"
        },
        "options":{
            "mom":""
        }
    }

    # Creates Z4 random wall at every tStep time slices starting at t0
    moduleTemplates["noiseRW"] = {
        "id":{
        "name":"",
            "type":"MSource::StagRandomWall"
        },
        "options":{
            "nSrc":"",  
            "tStep":"",
            "t0":"0",
            "colorDiag":"true",
        },
    }

    moduleTemplates["timeDilutedNoise"] = {
        "id":{
        "name":"",
            "type":"MNoise::StagTimeDilutedSpinColorDiagonal"
        },
        "options":{
            "nsrc":"",
            'tStep':'1'
        },
    }

    # Creates new vector of fields from given indices of chosen source vector
    moduleTemplates["splitVec"] = {
        "id":{
        "name":"",
            "type":"MUtilities::StagSourcePickIndices"
        },
        "options":{
            "source":"",  
            "indices":"",
        },
    }
    
    # Creates mixed precision CG solver
    moduleTemplates["cgMP"] = {
        "id":{
            "name":"",
            "type":"MSolver::StagMixedPrecisionCG"
        },
        "options":{
            "outerAction":"",
            "innerAction":"",
            "maxOuterIteration":"10000",
            "maxInnerIteration":"10000",
            "residual":"",
            "isEven":"false"
        },
    }

    # Projects chosen source onto low mode subspace of Dirac Matrix
    moduleTemplates["lmaProj"] = {
	"id":{
            "name":"",
            "type":"MSolver::StagLMA",
        },
        "options":{
            "action":"",
            "lowModes":"",
            "projector":"false",
            "eigStart":"0",
            "nEigs":"-1",
        }
    }

    # Applies chosen solver to sources
    moduleTemplates["quarkProp"] = {
        "id":{
            "name":"",
            "type":"MFermion::StagGaugeProp",
        },
        "options":{
            "source":"",
            "solver":"",
            "guess":"",
            "spinTaste":{
                "gammas":"",
                "gauge":""
            },
        },
    }

    # Contracts chosen quark propagators (q1, q2) to form correlator
    moduleTemplates["propContract"] = {
        "id":{
            "name":"",
            "type":"MContraction::StagMeson",
        },
        "options":{
            "source":"",
            "sink":"",
            "sinkFunc":"",
            "sourceShift":"",
            "sourceGammas":"",
            "sinkSpinTaste":{
                "gammas":"",
                "gauge":""
            },
            "output":"",
        },
    }

    # Builds A2A meson field
    moduleTemplates["mesonField"] = {
        "id":{
            "name":"",
            "type":"MContraction::StagA2AMesonField",
        },
        "options":{
            "action":"",
            "block":"",
            "mom":{
                "elem":"0 0 0",
            },
            "spinTaste":{
                "gammas":"",
                "gauge":""
            },
            "lowModes":"",
            "left":"",
            "right":"",
            "output":"",
        },
    }

    # Builds A2A Aslash meson field
    moduleTemplates["qedMesonField"] = {
        "id":{
            "name":"",
            "type":"MContraction::StagA2AASlashMesonField",
        },
        "options":{
            "action":"",
            "block":"",
            "mom":{
                "elem":"0 0 0",
            },
            "EmFunc":"",
            "nEmFields":"",
            "EmSeedString":"",
            "lowModes":"",
            "left":"",
            "right":"",
            "output":"",
        },
    }

    return moduleTemplates
