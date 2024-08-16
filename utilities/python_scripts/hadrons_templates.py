import copy

def generate_templates():
    module_templates = {}

    # Loads gauge fields from ILDG file
    module_templates["load_gauge"] = {
        "id":{
            "name":"",
            "type": "MIO::LoadIldg"
        },
        "options": {
            "file":""
        }
    }
    
    # Recasts double precision gauge field to single
    module_templates["cast_gauge"] = {
        "id":{
            "name":"",
            "type": "MUtilities::GaugeSinglePrecisionCast"
        },
        "options": {
            "field":""
        }
    }

    # Creates Dirac Matrix
    module_templates["action"] = {
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
    module_templates["action_float"] = copy.deepcopy(module_templates["action"])
    module_templates["action_float"]["id"]["type"] = "MAction::ImprovedStaggeredMILCF"

    module_templates["op"] = {
        "id":{
            "name":"",
            "type":"MFermion::StagOperators"
        },
        "options":{
            "action":""
        }
    }

    module_templates["irl"] = {
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
    module_templates["epack_load"] = {
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
    module_templates["eval_save"] = {
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
    module_templates["epack_modify"] = {
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
    module_templates["spin_taste"] = {
        "id":{
            "name":"",
            "type":"MFermion::SpinTaste"
        }
    }

    # Creates function that spacially ties up field with chosen momentum phase
    module_templates["sink"] = {
        "id":{
            "name":"",
            "type":"MSink::ScalarPoint"
        },
        "options":{
            "mom":""
        }
    }

    # Creates Z4 random wall at every time step time slices starting at t0
    module_templates["noise_rw"] = {
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

    module_templates["time_diluted_noise"] = {
        "id":{
        "name":"",
            "type":"MNoise::StagTimeDilutedSpinColorDiagonal"
        },
        "options":{
            "nsrc":"",
            'tStep':'1'
        },
    }

    module_templates["full_volume_noise"] = {
        "id":{
        "name":"",
            "type":"MNoise::StagFullVolumeSpinColorDiagonal"
        },
        "options":{
            "nsrc":""
        },
    }

    # Creates new vector of fields from given indices of chosen source vector
    module_templates["split_vec"] = {
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
    module_templates["mixed_precision_cg"] = {
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
    module_templates["lma_solver"] = {
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
    module_templates["quark_prop"] = {
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
    module_templates["prop_contract"] = {
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
    module_templates["meson_field"] = {
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
    module_templates["qed_meson_field"] = {
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

    # Create operator to generate EM gauge field (A_mu)
    module_templates["em_func"] = {
        "id":{
            "name":"",
            "type":"MGauge::StochEmFunc"
        },
        "options":{
            "gauge":"",
            "zmScheme":""
        }
    }

    # Create operator to generate EM gauge field (A_mu)
    module_templates["em_field"] = {
        "id":{
            "name":"",
            "type":"MGauge::StochEm"
        },
        "options":{
            "gauge":"",
            "zmScheme":"",
            "improvement":""
        }
    }

    module_templates["seq_aslash"] = {
        "id":{
        "name":"",
            "type":"MSource::StagSeqAslash"
        },
        "options":{
            "q":"",  
            "tA":"",
            "tB":"",
            "emField":"",
            "mom":"",
        },
    }

    module_templates["seq_gamma"] = {
        "id":{
        "name":"",
            "type":"MSource::StagSeqGamma"
        },
        "options":{
            "q":"",  
            "tA":"",
            "tB":"",
            "mom":"",
        },
    }

    # Builds A2A vectors
    module_templates["save_vector"] = {
        "id":{
            "name":"",
            "type":"MIO::SaveStagVector",
        },
        "options":{
            "field":"",
            "multiFile":"true",
            "output":""
        },
    }

    # Builds A2A vectors
    module_templates["a2a_vector"] = {
        "id":{
            "name":"",
            "type":"MSolver::StagA2AVectors",
        },
        "options":{
            "noise":"",
            "action":"",
            "lowModes":"",
            "solver":"",
            "highOutput":"",
            "norm2":"1.0",
            "highMultiFile":"false"
        },
    }
    
    # Loads A2A Vectors
    module_templates["load_vectors"] = {
        "id":{
            "name":"",
            "type":"MIO::StagLoadA2AVectors",
        },
        "options":{
            "filestem":"",
            "multiFile":"false",
            "size":""
        },
    }

    # A2A Contraction
    module_templates["contract_a2a_matrix"] = {
        "file":"",
        "dataset": "",
        "cacheSize":"1",
        "ni":"",
        "nj":"",
        "niOffset":"0",
        "njOffset":"0",
        "name":""
    }

    module_templates["contract_a2a_product"] = {
        "terms":"",
        "times": {
            "elem":"0"
        },
        "translations":"",
        "translationAverage":"true",
        "spaceNormalize":"false"
    }
        
    return module_templates
