#!/usr/bin/env python3
"""Consistency probe: stochastic A_slash A2A meson fields from the stencil
module (MContraction::StagA2AASlashMesonField, J.A local worker + CB pair
sources) vs the legacy module (MContraction::StagA2AASlashMesonFieldLegacy,
frozen A2AWorkerLocal + action/swapFn) on identical EmSeedString in the
same schedule (params/aslash-stencil-vs-legacy.20.xml, 4^4 exact IRL
eigenbasis).

The legacy path is NOT trusted ground truth (it has zero prior test
coverage): this script REPORTS the per-file max relative deviation.
deviation <= tol CORROBORATES both paths; larger deviation is an
investigation trigger for WHICH path is wrong -- never an automatic
stencil-side failure. Worker-level ground truth is Grid's
Test_a2a_stencil_ja_mom oracle suite.

Also structurally validates the stencil momentum stem (mf_stencil_mom):
file count nAmu*nMom, filename momentum == md.momentum metadata,
amu_index == filename channel token, dims equal the probe pair's.

Usage (from test/, after running the schedule
params/aslash-stencil-vs-legacy.20.xml with ../build-scalar/HadronsMILC
params/aslash-stencil-vs-legacy.20.xml --grid 4.4.4.4):

    python3 compare_aslash_mesonfield.py [--traj 20] [--tol 1e-6] [--strict]
--strict exits non-zero when any deviation exceeds tol (CI gate use).
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
STENCIL_STEM = "work/aslash-probe/m01/mf_stencil_cb"
LEGACY_STEM = "work/aslash-probe/m01/mf_legacy_cb"
MOM_STEM = "work/aslash-probe/m01/mf_stencil_mom"
N_EM = 2   # nEmFields
N_MOM = 2  # entries in the momentum stem's mom list


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


def max_rel_dev(a, b):
    den = np.maximum(np.abs(b), 1e-30)
    return float(np.max(np.abs(a - b) / den))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--traj", type=int, default=20)
    p.add_argument("--tol", type=float, default=1e-6,
                   help="corroboration tolerance (max relative deviation)")
    p.add_argument("--strict", action="store_true",
                   help="exit non-zero when deviation exceeds tol")
    args = p.parse_args()

    st_mats, st_mds = load(STENCIL_STEM, args.traj)
    lg_mats, lg_mds = load(LEGACY_STEM, args.traj)
    if set(st_mats) != set(lg_mats):
        sys.exit("file-name mismatch between stencil and legacy stems: "
                 "{} vs {}".format(sorted(st_mats), sorted(lg_mats)))

    worst, worst_name = 0.0, ""
    for name in sorted(st_mats):
        if st_mats[name].shape != lg_mats[name].shape:
            sys.exit("{}: shape mismatch {} vs {}".format(
                name, st_mats[name].shape, lg_mats[name].shape))
        dev = max_rel_dev(st_mats[name], lg_mats[name])
        print("  {:24s} max rel dev = {:.3e}".format(name, dev))
        if dev > worst:
            worst, worst_name = dev, name

    verdict = ("CORROBORATED (<= {}): the two implementations agree; both "
               "paths consistent".format(args.tol) if worst <= args.tol else
               "INVESTIGATE (> {}): stencil and legacy DISAGREE on {} -- "
               "determine WHICH path is wrong (legacy is not trusted "
               "ground truth; check Grid's Test_a2a_stencil_ja_mom first)"
               .format(args.tol, worst_name))
    print("consistency probe: max rel dev = {:.3e} -- {}".format(
        worst, verdict))
    if args.strict and worst > args.tol:
        sys.exit(1)

    # structural validation of the momentum stem
    mo_mats, mo_mds = load(MOM_STEM, args.traj)
    if len(mo_mats) != N_EM * N_MOM:
        sys.exit("momentum stem: expected {} files, got {}".format(
            N_EM * N_MOM, len(mo_mats)))
    for name, md in mo_mds.items():
        tok = name.split("_")  # seed_a_p1_p2_p3
        if int(tok[-4]) != md["amu_index"]:
            sys.exit("{}: filename channel {} != amu_index {}".format(
                name, tok[-4], md["amu_index"]))
        fmom = [float(t) for t in tok[-3:]]
        if not np.array_equal(np.asarray(fmom), md["momentum"]):
            sys.exit("{}: filename momentum {} != metadata {}".format(
                name, fmom, md["momentum"].tolist()))
    ref_shape = next(iter(st_mats.values())).shape
    for name, m in mo_mats.items():
        if m.shape != ref_shape:
            sys.exit("{}: shape {} != {}".format(name, m.shape, ref_shape))
    print("momentum stem structure OK: {} files, names/metadata/dims "
          "consistent".format(len(mo_mats)))


if __name__ == "__main__":
    main()
