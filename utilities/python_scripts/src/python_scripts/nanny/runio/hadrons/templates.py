import typing as t

def xml_wrapper(runid: str, sched: str, cfg: str) -> dict:
    params = {
        "grid": {
            "parameters": {
                "runId": runid,
                "trajCounter": {
                    "start": cfg,
                    "end": "10000",
                    "step": "10000",
                },
                "genetic": {
                    "popSize": "20",
                    "maxGen": "1000",
                    "maxCstGen": "100",
                    "mutationRate": "0.1",
                },
                "graphFile": "",
                "scheduleFile": sched,
                "saveSchedule": "false",
                "parallelWriteMaxRetry": "-1",
            },
            "modules": {},
        },
    }

    return params


def load_gauge(name: str, file: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MIO::LoadIldg"
        },
        "options": {
            "file": file
        }
    }


def cast_gauge(name: str, field: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MUtilities::GaugeSinglePrecisionCast"
        },
        "options": {
            "field": field
        }
    }


def action(name: str, mass: str, gauge_fat: str, gauge_long: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MAction::ImprovedStaggeredMILC"
        },
        "options": {
            "mass": mass,
            "c1": "1.0",
            "c2": "1.0",
            "tad": "1.0",
            "boundary": "1 1 1 1",
            "twist": "0 0 0",
            "Ls": "1",
            "gaugefat": gauge_fat,
            "gaugelong": gauge_long
        }
    }


def action_float(*args, **kwargs) -> t.Dict:
    res = action(*args, **kwargs)
    res['id']['type'] = "MAction::ImprovedStaggeredMILCF"
    return res


def op(name: str, action: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MFermion::StagOperators"
        },
        "options": {
            "action": action
        }
    }


def irl(name: str,
        op: str,
        alpha: str,
        beta: str,
        npoly: str,
        nstop: str,
        nk: str,
        nm: str,
        multifile: str,
        output: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MSolver::StagFermionIRL"
        },
        "options": {
            "op": op,
            "lanczosParams": {
                "Cheby": {
                    "alpha": alpha,
                    "beta": beta,
                    "Npoly": npoly

                },
                "Nstop": nstop,
                "Nk": nk,
                "Nm": nm,
                "resid": "1e-8",
                "MaxIt": "5000",
                "betastp": "0",
                "MinRes": "0"
            },
            "evenEigen": "false",
            "redBlack": "true",
            "output": output,
            "multiFile": multifile
        }
    }


def epack_load(name: str, filestem: str, size: str, multifile: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MIO::StagLoadFermionEigenPack"
        },
        "options": {
            "redBlack": "true",
            "filestem": filestem,
            "multiFile": multifile,
            "size": size,
            "Ls": "1"
        }
    }


def eval_save(name: str, eigen_pack: str, output: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MUtilities::EigenPackExtractEvals"
        },
        "options": {
            "eigenPack": eigen_pack,
            "output": output

        }
    }


def epack_modify(name: str, eigen_pack: str, mass: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MUtilities::ModifyEigenPackMILC"
        },
        "options": {
            "eigenPack": eigen_pack,
            "mass": mass,
            "evenEigen": "false",
            "normalizeCheckerboard": "false"
        }
    }


def spin_taste(name: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MFermion::SpinTaste"
        }
    }


def sink(name: str, mom: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MSink::ScalarPoint"
        },
        "options": {
            "mom": mom

        }
    }


def noise_rw(name: str, nsrc: str, t0: str, tstep: str, noise: str='') -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MSource::StagRandomWall"
        },
        "options": {
            "nSrc": nsrc,
            "tStep": tstep,
            "t0": t0,
            "colorDiag": "true",
            'noise':noise
        },
    }


def time_diluted_noise(name: str, nsrc: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MNoise::StagTimeDilutedSpinColorDiagonal"
        },
        "options": {
            "nsrc": nsrc,
            'tStep': '1'
        },
    }


def full_volume_noise(name: str, nsrc: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MNoise::StagFullVolumeSpinColorDiagonal"
        },
        "options": {
            "nsrc": nsrc

        },
    }


def split_vec(name: str, source: str, indices: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MUtilities::StagSourcePickIndices"
        },
        "options": {
            "source": source,
            "indices": indices,
        },
    }


def mixed_precision_cg(name: str, outer_action: str, inner_action: str, residual: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MSolver::StagMixedPrecisionCG"
        },
        "options": {
            "outerAction": outer_action,
            "innerAction": inner_action,
            "maxOuterIteration": "10000",
            "maxInnerIteration": "10000",
            "residual": residual,
            "isEven": "false"
        },
    }


def lma_solver(name: str, action: str, low_modes: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MSolver::StagLMA",
        },
        "options": {
            "action": action,
            "lowModes": low_modes,
            "projector": "false",
            "eigStart": "0",
            "nEigs": "-1",
        }
    }


def quark_prop(name: str,
               source: str,
               solver: str,
               guess: str,
               gammas: str,
               gauge: str,
               apply_g5: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MFermion::StagGaugeProp",
        },
        "options": {
            "source": source,
            "solver": solver,
            "guess": guess,
            "spinTaste": {
                "gammas": gammas,
                "gauge": gauge,
                "applyG5":apply_g5

            },
        },
    }


def prop_contract(name: str,
                  source: str,
                  sink: str,
                  sink_func: str,
                  source_shift: str,
                  source_gammas: str,
                  sink_gammas: str,
                  gauge: str,
                  output: str,
                  apply_g5: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MContraction::StagMeson",
        },
        "options": {
            "source": source,
            "sink": sink,
            "sinkFunc": sink_func,
            "sourceShift": source_shift,
            "sourceGammas": source_gammas,
            "sinkSpinTaste": {
                "gammas": sink_gammas,
                "gauge": gauge,
                "applyG5": apply_g5
            },
            "output": output,
        },
    }


def meson_field(name: str,
                action: str,
                block: str,
                gammas: str,
                gauge: str,
                low_modes: str,
                left: str,
                right: str,
                output: str,
                apply_g5:str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MContraction::StagA2AMesonField",
        },
        "options": {
            "action": action,
            "block": block,
            "mom": {
                "elem": "0 0 0",
            },
            "spinTaste": {
                "gammas": gammas,
                "gauge": gauge,
                "applyG5": apply_g5

            },
            "lowModes": low_modes,
            "left": left,
            "right": right,
            "output": output,
        },
    }


def qed_meson_field(name: str,
                    action: str,
                    block: str,
                    em_func: str,
                    n_em_fields: str,
                    em_seed_string: str,
                    low_modes: str,
                    left: str,
                    right: str,
                    output: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MContraction::StagA2AASlashMesonField",
        },
        "options": {
            "action": action,
            "block": block,
            "mom": {
                "elem": "0 0 0",
            },
            "EmFunc": em_func,
            "nEmFields": n_em_fields,
            "EmSeedString": em_seed_string,
            "lowModes": low_modes,
            "left": left,
            "right": right,
            "output": output,
        },
    }


def em_func(name: str, gauge: str, zm_scheme: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MGauge::StochEmFunc"
        },
        "options": {
            "gauge": gauge,
            "zmScheme": zm_scheme

        }
    }


def em_field(name: str, gauge: str, zm_scheme: str, improvement: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MGauge::StochEm"
        },
        "options": {
            "gauge": gauge,
            "zmScheme": zm_scheme,
            "improvement": improvement

        }
    }


def seq_aslash(name: str, q: str, ta: str, tb: str, em_field: str, mom: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MSource::StagSeqAslash"
        },
        "options": {
            "q": q,
            "tA": ta,
            "tB": tb,
            "emField": em_field,
            "mom": mom,
        },
    }


def seq_gamma(name: str, q: str, ta: str, tb: str, mom: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MSource::StagSeqGamma"
        },
        "options": {
            "q": q,
            "tA": ta,
            "tB": tb,
            "mom": mom,
        },
    }


def save_vector(name: str, field: str, output: str, multifile: str='false') -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MIO::SaveStagVector",
        },
        "options": {
            "field": field,
            "multiFile": multifile,
            "output": output

        },
    }


def a2a_vector(name: str, noise: str, action: str, low_modes: str, solver: str, high_output: str) -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MSolver::StagA2AVectors",
        },
        "options": {
            "noise": noise,
            "action": action,
            "lowModes": low_modes,
            "solver": solver,
            "highOutput": high_output,
            "norm2": "1.0",
            "highMultiFile": "false"
        },
    }


def load_vectors(name: str, filestem: str, size: str, multifile: str = 'false') -> t.Dict:
    return {
        "id": {
            "name": name,
            "type": "MIO::StagLoadA2AVectors",
        },
        "options": {
            "filestem": filestem,
            "multiFile": multifile,
            "size": size

        },
    }