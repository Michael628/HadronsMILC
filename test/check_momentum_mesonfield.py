#!/usr/bin/env python3
"""Validate momentum-projected A2A meson fields from
MContraction::StagA2AMesonField against the zero-momentum control module in
the same schedule (params/a2a-momentum-irl-exact.20.xml, 4^4 exact IRL
eigenbasis).

Reference-free wiring gate for the g-major flat-mu dims swap:
  1. file counts: control = 1 file per gamma; momentum module = nMom files
     per gamma (GammaName_p1_p2_p3.h5)
  2. the momentum module's zero-momentum file is VALUE-IDENTICAL to the
     control module's file for each gamma (exp(0)=1 exactly: the p=0
     projection reduces to the unprojected sum with an identical
     accumulation order -- a transposed (momentum, gamma) allocation would
     place wrong-pair content in every file and fail this check)
  3. per file: momentum parsed from the filename == md.momentum metadata;
     the gamma metadata pair (gamma_spin/gamma_taste attrs) matches the
     control file carrying the same gamma-name prefix; dataset dims equal
     the control's

Usage (from test/, after running the schedule
params/a2a-momentum-irl-exact.20.xml with ../build-scalar/HadronsMILC
params/a2a-momentum-irl-exact.20.xml --grid 4.4.4.4):

    python3 check_momentum_mesonfield.py [--traj 20] [--tol 0.0]
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
CTL_STEM = "work/e384n1dt1-momtest/mesons/m01/mf_ctl"
MOM_STEM = "work/e384n1dt1-momtest/mesons/m01/mf_mom"
N_MOM = 6  # entries in the momentum module's mom list
ND_MOM = 3  # momentum components per entry (nd-1)


def load(stem, traj):
    outdir = "{}.{}".format(os.path.join(TEST_DIR, stem), traj)
    files = sorted(glob.glob(os.path.join(outdir, "*.h5")))
    if not files:
        sys.exit("no HDF5 files under {}".format(outdir))
    mats, mds = {}, {}
    for path in files:
        name = os.path.basename(path)[:-3]
        with h5py.File(path, "r") as f:
            d = f[name]["a2aMatrix"][()]
            mats[name] = (d["re"] + 1j * d["im"]).astype(np.complex64)
            grp = f[name]
            mom = np.asarray(grp.attrs["momentum"], dtype=float)
            spin = np.asarray(grp.attrs["gamma_spin"])
            taste = np.asarray(grp.attrs["gamma_taste"])
            mds[name] = (mom, spin, taste)
    return mats, mds


def split_name(name):
    """'G1_G1_1_0_0' -> ('G1_G1', [1.0, 0.0, 0.0]) (last ND_MOM tokens are
    the momentum components)."""
    tok = name.split("_")
    mom = [float(t) for t in tok[-ND_MOM:]]
    gamma = "_".join(tok[:-ND_MOM])
    return gamma, mom


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--traj", type=int, default=20)
    p.add_argument("--tol", type=float, default=0.0,
                   help="max allowed |delta| for the p=0 identity "
                        "(default 0: bitwise; a small nonzero value "
                        "indicates rounding-order drift to investigate)")
    args = p.parse_args()

    ctl_mats, ctl_mds = load(CTL_STEM, args.traj)
    mom_mats, mom_mds = load(MOM_STEM, args.traj)

    # 1. file counts
    n_gamma = len(ctl_mats)
    if len(mom_mats) != n_gamma * N_MOM:
        sys.exit("expected {} momentum files ({} gammas x {} momenta), "
                 "got {}".format(n_gamma * N_MOM, n_gamma, N_MOM,
                                 len(mom_mats)))
    print("file counts OK: {} control, {} momentum".format(
        n_gamma, len(mom_mats)))

    ctl_by_gamma = {}
    for name, (mom, spin, taste) in ctl_mds.items():
        gamma, fmom = split_name(name)
        ctl_by_gamma[gamma] = (name, mom, spin, taste)

    # 2 + 3. per-file checks
    worst = 0.0
    for name, (mom, spin, taste) in mom_mds.items():
        gamma, fmom = split_name(name)
        # filename momentum == metadata momentum
        if not np.array_equal(np.asarray(fmom), mom):
            sys.exit("{}: filename momentum {} != metadata {}".format(
                name, fmom, mom.tolist()))
        # gamma metadata matches the control file with the same prefix
        ctl_name, ctl_mom, ctl_spin, ctl_taste = ctl_by_gamma[gamma]
        if (spin, taste) != (ctl_spin, ctl_taste):
            sys.exit("{}: gamma metadata ({}, {}) != control ({}, "
                     "{})".format(name, spin, taste, ctl_spin, ctl_taste))
        # dims equal the control's
        if mom_mats[name].shape != ctl_mats[ctl_name].shape:
            sys.exit("{}: shape {} != control shape {}".format(
                name, mom_mats[name].shape, ctl_mats[ctl_name].shape))
        # the zero-momentum files are value-identical to the control
        if np.all(mom == 0.0):
            dev = np.max(np.abs(mom_mats[name] - ctl_mats[ctl_name]))
            worst = max(worst, dev)
            if dev > args.tol:
                sys.exit("{}: p=0 identity violated, max |delta| = {} "
                         "(> tol {})".format(name, dev, args.tol))
    print("momentum wiring OK: filenames, metadata, dims consistent; "
          "p=0 identity max |delta| = {}".format(worst))


if __name__ == "__main__":
    main()
