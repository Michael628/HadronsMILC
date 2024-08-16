#! /usr/bin/env python3

# Python 3 version

import sys, os, yaml, re, subprocess, copy
from todo_utils import *

def main():
    params = load_param("params.yaml")["LMIparam"]

    ens, eigs, noise, dt, start, stop, mass = (params[k] for k in ["ENS","EIGS","NOISE","DT","TSTART","TSTOP","TESTMASSLMA"])

    massa2a=params["TESTMASSA2A"]

    configs = sys.argv[1:]

    if configs[0] == "a2a":
        for series,cfg in [c.split(".") for c in configs[1:]]:
            for i,mass1 in enumerate(params["MASSES"].split()):
                for j,mass2 in enumerate(params["MASSES"].split()):
                    if i == j:
                        m = "m"+mass1
                    elif j > i:
                        continue
                    else:
                        m = "m"+mass1+"_m"+mass2

                    for k,gamma in enumerate(["pion","vecX","vecY","vecZ","vecX_onelink","vecY_onelink","vecZ_onelink"]):
                        if i != 0 and k > 3:
                            continue
                        if i != j:
                            filename = f"e{eigs}n{noise}dt{dt}/correlators/{m}/a2aLL/series_{series}/{gamma}1_0_{gamma}2.{cfg}.h5"
                        else:
                            filename = f"e{eigs}n{noise}dt{dt}/correlators/{m}/a2aLL/series_{series}/{gamma}_0_{gamma}.{cfg}.h5"
                        if not os.path.isfile(filename):
                            print(f"Missing file: {filename}")
                        
        exit(0)
    for series,cfg in [c.split(".") for c in configs]:
        for lat in ["l","lng","fat"]:
            filename = f"configs/{lat}{ens}{series}.ildg.{cfg}"
            if not os.path.isfile(filename):
                print(f"Missing file: {filename}")
                              
        for gamma in ["pion_local","vec_local","vec_onelink"]:
            m= mass if "local" in gamma else massa2a
            for dset in ["ama","ranLL"]:
                for t in range(int(start),int(stop)+1,int(dt)):
                    filename = f"e{eigs}n{noise}dt{dt}/correlators/{gamma}/{dset}/corr_{gamma}_{dset}_m{m}_t{str(t)}_{series}.{cfg}.h5"
                    if not os.path.isfile(filename):
                        filename = f"e{eigs}n{noise}dt{dt}/correlators/m{m}/{gamma}/{dset}/corr_{gamma}_{dset}_m{m}_t{str(t)}_{series}.{cfg}.h5"
                        if not os.path.isfile(filename):
                            print(f"Missing file: {filename}")
                        
        for gamma in ["G5_G5","GX_GX","GY_GY","GZ_GZ","GX_G1","GY_G1","GZ_G1"]:
            filename = f"e{eigs}n{noise}dt{dt}/mesons/mf_{series}.{cfg}/{gamma}_0_0_0.h5"
            if not os.path.isfile(filename):
                filename = f"e{eigs}n{noise}dt{dt}/mesons/m{massa2a}/mf_{series}.{cfg}/{gamma}_0_0_0.h5"
                if not os.path.isfile(filename):
                    print(f"Missing file: {filename}")

if __name__ == "__main__":
    main()

