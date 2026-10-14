#!/usr/bin/env python3
"""Structural validation of one-link J.A A2A meson fields
(MContraction::StagA2AASlashOneLinkMesonField, params/aslash-onelink.20.xml,
4^4 exact IRL eigenbasis): file count nSample*4*nMom; per file the
filename momentum == md.momentum, filename component token == amu_index,
md.amu_seed == "<seed>_<sample>"; dims consistent; and non-degeneracy
(different samples/components/momenta produce different data -- catches a
stuck generation loop or a wrong sample/channel slot mapping).

Numeric ground truth for the one-link operator itself is Grid's
Test_a2a_stencil_ja_mom oracle suite (worker-level); this script gates the
HadronsMILC-side wiring only.

Usage (from test/, after running the schedule with
../build-scalar/HadronsMILC params/aslash-onelink.20.xml --grid 4.4.4.4):

    python3 check_aslash_onelink.py [--traj 20]
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
STEM = "work/ja1link/m01/mf"
N_SAMPLE = 2  # nEmFields
N_COMP = 4    # J_kA_k channels (Nd)
N_MOM = 2     # mom entries
SEED = "ja1link-seed"


def load(traj):
    outdir = "{}.{}".format(os.path.join(TEST_DIR, STEM), traj)
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
            seed = np.asarray(grp.attrs["amu_seed"]).ravel()[0]
            if isinstance(seed, bytes):
                seed = seed.decode()
            mds[name] = {
                "momentum": np.asarray(grp.attrs["momentum"], dtype=float),
                "amu_seed": str(seed),
                "amu_index": int(
                    np.asarray(grp.attrs["amu_index"]).ravel()[0]),
            }
    return mats, mds


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--traj", type=int, default=20)
    args = p.parse_args()

    mats, mds = load(args.traj)
    expect = N_SAMPLE * N_COMP * N_MOM
    if len(mats) != expect:
        sys.exit("expected {} files ({}x{}x{}), got {}".format(
            expect, N_SAMPLE, N_COMP, N_MOM, len(mats)))

    ref_shape = None
    for name, md in mds.items():
        tok = name.split("_")  # seed_j_k_p1_p2_p3
        j, k = int(tok[-5]), int(tok[-4])
        fmom = [float(t) for t in tok[-3:]]
        if k != md["amu_index"]:
            sys.exit("{}: filename component {} != amu_index {}".format(
                name, k, md["amu_index"]))
        if md["amu_seed"] != "{}_{}".format(SEED, j):
            sys.exit("{}: amu_seed {} != expected {}_{}".format(
                name, md["amu_seed"], SEED, j))
        if not np.array_equal(np.asarray(fmom), md["momentum"]):
            sys.exit("{}: filename momentum {} != metadata {}".format(
                name, fmom, md["momentum"].tolist()))
        if ref_shape is None:
            ref_shape = mats[name].shape
        elif mats[name].shape != ref_shape:
            sys.exit("{}: shape {} != {}".format(
                name, mats[name].shape, ref_shape))

    # non-degeneracy: all files pairwise distinct
    names = sorted(mats)
    for i in range(len(names)):
        for l in range(i + 1, len(names)):
            if np.array_equal(mats[names[i]], mats[names[l]]):
                sys.exit("degenerate pair: {} == {} bitwise".format(
                    names[i], names[l]))

    print("one-link structure OK: {} files, naming/metadata/dims "
          "consistent, all distinct".format(len(mats)))


if __name__ == "__main__":
    main()
