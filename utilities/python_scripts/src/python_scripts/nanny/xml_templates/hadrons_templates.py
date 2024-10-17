from python_scripts.utils import ReadOnlyDict


module_templates = {}

# Loads gauge fields from ILDG file
module_templates["load_gauge"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MIO::LoadIldg"
    }),
    "options": ReadOnlyDict({
        "file": ""
    })
})

# Recasts double precision gauge field to single
module_templates["cast_gauge"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MUtilities::GaugeSinglePrecisionCast"
    }),
    "options": ReadOnlyDict({
        "field": ""
    })
})

# Creates Dirac Matrix
module_templates["action"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MAction::ImprovedStaggeredMILC"
    }),
    "options": ReadOnlyDict({
        "mass": "",
        "c1": "1.0",
        "c2": "1.0",
        "tad": "1.0",
        "boundary": "1 1 1 1",
        "twist": "0 0 0",
        "Ls": "1",
        "gaugefat": "",
        "gaugelong": ""
    })
})

# Creates single precision Dirac Matrix
module_templates["action_float"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MAction::ImprovedStaggeredMILCF"
    }),
    "options": module_templates['action']['options']
})

module_templates["op"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MFermion::StagOperators"
    }),
    "options": ReadOnlyDict({
        "action": ""
    })
})

module_templates["irl"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MSolver::StagFermionIRL"
    }),
    "options": ReadOnlyDict({
        "op": "",
        "lanczosParams": ReadOnlyDict({
            "Cheby": ReadOnlyDict({
                "alpha": "",
                "beta": "",
                "Npoly": ""
            }),
            "Nstop": "",
            "Nk": "",
            "Nm": "",
            "resid": "1e-8",
            "MaxIt": "5000",
            "betastp": "0",
            "MinRes": "0"
        }),
        "evenEigen": "false",
        "redBlack": "true",
        "output": "",
        "multiFile": "false"
    })
})

# Loads eigenvector pack
module_templates["epack_load"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MIO::StagLoadFermionEigenPack"
    }),
    "options": ReadOnlyDict({
        "redBlack": "true",
        "filestem": "",
        "multiFile": "false",
        "size": "",
        "Ls": "1"
    })
})

# Write eigenvalues into separate XML file
module_templates["eval_save"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MUtilities::EigenPackExtractEvals"
    }),
    "options": ReadOnlyDict({
        "eigenPack": "",
        "output": ""
    })
})

# Shift eigenvalues of eigenpack by a mass value
module_templates["epack_modify"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MUtilities::ModifyEigenPackMILC"
    }),
    "options": ReadOnlyDict({
        "eigenPack": "",
        "mass": "",
        "evenEigen": "false",
        "normalizeCheckerboard": "false"
    })
})

# Create operator to apply spin taste to fields
module_templates["spin_taste"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MFermion::SpinTaste"
    })
})

# Creates function that spacially ties up field with chosen momentum phase
module_templates["sink"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MSink::ScalarPoint"
    }),
    "options": ReadOnlyDict({
        "mom": ""
    })
})

# Creates Z4 random wall at every time step time slices starting at t0
module_templates["noise_rw"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MSource::StagRandomWall"
    }),
    "options": ReadOnlyDict({
        "nSrc": "",
        "tStep": "",
        "t0": "0",
        "colorDiag": "true",
    }),
})

module_templates["time_diluted_noise"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MNoise::StagTimeDilutedSpinColorDiagonal"
    }),
    "options": ReadOnlyDict({
        "nsrc": "",
        'tStep': '1'
    }),
})

module_templates["full_volume_noise"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MNoise::StagFullVolumeSpinColorDiagonal"
    }),
    "options": ReadOnlyDict({
        "nsrc": ""
    }),
})

# Creates new vector of fields from given indices of chosen source vector
module_templates["split_vec"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MUtilities::StagSourcePickIndices"
    }),
    "options": ReadOnlyDict({
        "source": "",
        "indices": "",
    }),
})

# Creates mixed precision CG solver
module_templates["mixed_precision_cg"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MSolver::StagMixedPrecisionCG"
    }),
    "options": ReadOnlyDict({
        "outerAction": "",
        "innerAction": "",
        "maxOuterIteration": "10000",
        "maxInnerIteration": "10000",
        "residual": "",
        "isEven": "false"
    }),
})

# Projects chosen source onto low mode subspace of Dirac Matrix
module_templates["lma_solver"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MSolver::StagLMA",
    }),
    "options": ReadOnlyDict({
        "action": "",
        "lowModes": "",
        "projector": "false",
        "eigStart": "0",
        "nEigs": "-1",
    })
})

# Applies chosen solver to sources
module_templates["quark_prop"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MFermion::StagGaugeProp",
    }),
    "options": ReadOnlyDict({
        "source": "",
        "solver": "",
        "guess": "",
        "spinTaste": ReadOnlyDict({
            "gammas": "",
            "gauge": ""
        }),
    }),
})

# Contracts chosen quark propagators (q1, q2) to form correlator
module_templates["prop_contract"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MContraction::StagMeson",
    }),
    "options": ReadOnlyDict({
        "source": "",
        "sink": "",
        "sinkFunc": "",
        "sourceShift": "",
        "sourceGammas": "",
        "sinkSpinTaste": ReadOnlyDict({
            "gammas": "",
            "gauge": ""
        }),
        "output": "",
    }),
})

# Builds A2A meson field
module_templates["meson_field"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MContraction::StagA2AMesonField",
    }),
    "options": ReadOnlyDict({
        "action": "",
        "block": "",
        "mom": ReadOnlyDict({
            "elem": "0 0 0",
        }),
        "spinTaste": ReadOnlyDict({
            "gammas": "",
            "gauge": ""
        }),
        "lowModes": "",
        "left": "",
        "right": "",
        "output": "",
    }),
})

# Builds A2A Aslash meson field
module_templates["qed_meson_field"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MContraction::StagA2AASlashMesonField",
    }),
    "options": ReadOnlyDict({
        "action": "",
        "block": "",
        "mom": ReadOnlyDict({
            "elem": "0 0 0",
        }),
        "EmFunc": "",
        "nEmFields": "",
        "EmSeedString": "",
        "lowModes": "",
        "left": "",
        "right": "",
        "output": "",
    }),
})

# Create operator to generate EM gauge field (A_mu)
module_templates["em_func"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MGauge::StochEmFunc"
    }),
    "options": ReadOnlyDict({
        "gauge": "",
        "zmScheme": ""
    })
})

# Create operator to generate EM gauge field (A_mu)
module_templates["em_field"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MGauge::StochEm"
    }),
    "options": ReadOnlyDict({
        "gauge": "",
        "zmScheme": "",
        "improvement": ""
    })
})

module_templates["seq_aslash"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MSource::StagSeqAslash"
    }),
    "options": ReadOnlyDict({
        "q": "",
        "tA": "",
        "tB": "",
        "emField": "",
        "mom": "",
    }),
})

module_templates["seq_gamma"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MSource::StagSeqGamma"
    }),
    "options": ReadOnlyDict({
        "q": "",
        "tA": "",
        "tB": "",
        "mom": "",
    }),
})


# Builds A2A vectors
module_templates["save_vector"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MIO::SaveStagVector",
    }),
    "options": ReadOnlyDict({
        "field": "",
        "multiFile": "true",
        "output": ""
    }),
})

# Builds A2A vectors
module_templates["a2a_vector"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MSolver::StagA2AVectors",
    }),
    "options": ReadOnlyDict({
        "noise": "",
        "action": "",
        "lowModes": "",
        "solver": "",
        "highOutput": "",
        "norm2": "1.0",
        "highMultiFile": "false"
    }),
})

# Loads A2A Vectors
module_templates["load_vectors"] = ReadOnlyDict({
    "id": ReadOnlyDict({
        "name": "",
        "type": "MIO::StagLoadA2AVectors",
    }),
    "options": ReadOnlyDict({
        "filestem": "",
        "multiFile": "false",
        "size": ""
    }),
})

# A2A Contraction
module_templates["contract_a2a_matrix"] = ReadOnlyDict({
    "file": "",
    "dataset": "",
    "cacheSize": "1",
    "ni": "",
    "nj": "",
    "niOffset": "0",
    "njOffset": "0",
    "name": ""
})

module_templates["contract_a2a_product"] = ReadOnlyDict({
    "terms": "",
    "times": ReadOnlyDict({
        "elem": "0"
    }),
    "translations": "",
    "translationAverage": "true",
    "spaceNormalize": "false"
})

module_templates = ReadOnlyDict(module_templates)
