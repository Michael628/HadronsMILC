#! /usr/bin/env python3

# Python 3 version

import sys
import os
import python_scripts.utils as utils


def main():
    params = utils.load_param("params.yaml")["lmi_param"]

    ens, eigs, noise, dt, start, stop, mass = (
        params[k] for k in ["ENS", "EIGS", "NOISE",
                            "DT", "TSTART", "TSTOP", "TESTMASSLMA"])

    massa2a = params["TESTMASSA2A"]

    configs = sys.argv[1:]

    if configs[0] == "a2a":
        for series, cfg in [c.split(".") for c in configs[1:]]:
            for i, mass1 in enumerate(params["MASSES"].split()):
                for j, mass2 in enumerate(params["MASSES"].split()):
                    if i == j:
                        m = "m" + mass1
                    elif j > i:
                        continue
                    else:
                        m = "m" + mass1 + "_m" + mass2

                    for k, gamma in enumerate(
                        ["pion",
                         "vecX",
                         "vecY",
                         "vecZ",
                         "vecX_onelink",
                         "vecY_onelink",
                         "vecZ_onelink"]):
                        if i != 0 and k > 3:
                            continue
                        if i != j:
                            filename = (
                                f"e{eigs}n{noise}dt{dt}/correlators/"
                                f"{m}/a2aLL/series_{series}/"
                                f"{gamma}1_0_{gamma}2.{cfg}.h5"
                            )
                        else:
                            filename = (
                                f"e{eigs}n{noise}dt{dt}/correlators/"
                                f"{m}/a2aLL/series_{series}/"
                                f"{gamma}_0_{gamma}.{cfg}.h5"
                            )
                        if not os.path.isfile(filename):
                            print(f"Missing file: {filename}")

        exit(0)
    for series, cfg in [c.split(".") for c in configs]:
        for lat in ["l", "lng", "fat"]:
            filename = f"configs/{lat}{ens}{series}.ildg.{cfg}"
            if not os.path.isfile(filename):
                print(f"Missing file: {filename}")

        for gamma in ["pion_local", "vec_local", "vec_onelink"]:
            m = mass if "local" in gamma else massa2a
            for dset in ["ama", "ranLL"]:
                for t in range(int(start), int(stop) + 1, int(dt)):
                    filename = (
                        f"e{eigs}n{noise}dt{dt}/correlators/"
                        f"{gamma}/{dset}/"
                        f"corr_{gamma}_{dset}_m{m}_t{str(t)}_{series}.{cfg}.h5"
                    )
                    if not os.path.isfile(filename):
                        filename = (
                            f"e{eigs}n{noise}dt{dt}/correlators/"
                            f"m{m}/{gamma}/{dset}/corr_{gamma}_{dset}"
                            f"_m{m}_t{str(t)}_{series}.{cfg}.h5"
                        )
                        if not os.path.isfile(filename):
                            print(f"Missing file: {filename}")

        for gamma in ["G5_G5", "GX_GX", "GY_GY",
                      "GZ_GZ", "GX_G1", "GY_G1", "GZ_G1"]:
            filename = (
                f"e{eigs}n{noise}dt{dt}/mesons/"
                f"mf_{series}.{cfg}/{gamma}_0_0_0.h5"
            )
            if not os.path.isfile(filename):
                filename = (
                    f"e{eigs}n{noise}dt{dt}/mesons/m{massa2a}/"
                    f"mf_{series}.{cfg}/{gamma}_0_0_0.h5"
                )
                if not os.path.isfile(filename):
                    print(f"Missing file: {filename}")


if __name__ == "__main__":
    main()
