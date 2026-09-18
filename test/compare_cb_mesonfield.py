#!/usr/bin/env python3
"""Numerically compare checkerboarded-low-mode meson fields written by the
stencil module (MContraction::StagA2AMesonField with cbPairsLeft/cbPairsRight)
against the legacy module (MContraction::StagA2AMesonFieldLegacy with action)
on the same exact IRL eigenbasis (4^4 free-field configuration, 384 odd-CB
eigenvectors = the full preconditioned space).

Both modules consume the identical mass-shifted checkerboarded eigenpack; the
stencil path packs CB pairs via MUtilities::EigenPackCBPairs instances and
converts the worker's raw parity partials to the legacy interleaved (M, M-dagger)
slot layout, so every a2aMatrix dataset must agree to float32 rounding. The
schedule writes popcount-SEPARATED module pairs (the legacy module routes all
popcount>=2 gammas through one A2AWorkerSpinTaste, which Grid constrains to a
uniform popcount), spanning popcount 0-4 and both sigma signs of the
reconstruction table. Precedent: compare_lma_mesonfield.py.

Usage (from the test/ directory, after running the schedule
params/a2a-cb-mesonfields-irl-exact.20.xml with ../HadronsMILC --grid 4.4.4.4):

    python3 compare_cb_mesonfield.py [--traj 20] [--tol 1e-6]

Pass --stencil/--legacy (together) to compare a single module pair instead of
all default pairs.
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
# (stencil stem, legacy stem) per popcount-separated module pair, matching
# the module names/output stems in params/a2a-cb-mesonfields-irl-exact.20.xml
STEM_PAIRS = [
    ("e384n1dt1-cb/mesons/m01/mf_stencil_pc012",
     "e384n1dt1-cb/mesons/m01/mf_legacy_pc012"),
    ("e384n1dt1-cb/mesons/m01/mf_stencil_pc3",
     "e384n1dt1-cb/mesons/m01/mf_legacy_pc3"),
    ("e384n1dt1-cb/mesons/m01/mf_stencil_pc4",
     "e384n1dt1-cb/mesons/m01/mf_legacy_pc4"),
]


def load(stem, traj):
    outdir = "{}.{}".format(os.path.join(TEST_DIR, stem), traj)
    files = sorted(glob.glob(os.path.join(outdir, "*.h5")))
    if not files:
        sys.exit("no HDF5 files under {}".format(outdir))
    mats = {}
    for path in files:
        name = os.path.basename(path)[:-3]
        with h5py.File(path, "r") as f:
            d = f[name]["a2aMatrix"][()]
        # the writer stores an HDF5 compound type (re/im float pair), not a
        # native complex type -- view it as complex64 before any arithmetic
        mats[name] = (d["re"] + 1j * d["im"]).astype(np.complex64)
    return mats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stencil", default=None,
                   help="single stencil stem override (needs --legacy too)")
    p.add_argument("--legacy", default=None,
                   help="single legacy stem override (needs --stencil too)")
    p.add_argument("--traj", type=int, default=20)
    p.add_argument("--tol", type=float, default=1e-6,
                   help="max relative deviation (default 1e-6: the two "
                        "kernels differ only in FP rounding order, output "
                        "is float32)")
    args = p.parse_args()

    if (args.stencil is None) != (args.legacy is None):
        sys.exit("--stencil and --legacy must be given together")
    pairs = ([(args.stencil, args.legacy)] if args.stencil is not None
             else STEM_PAIRS)

    worst = 0.0
    for stencil_stem, legacy_stem in pairs:
        stn = load(stencil_stem, args.traj)
        leg = load(legacy_stem, args.traj)

        if set(stn) != set(leg):
            sys.exit("dataset name mismatch: stencil {} vs legacy {}".format(
                sorted(stn), sorted(leg)))

        print("== {} vs {} ==".format(
            os.path.basename(stencil_stem), os.path.basename(legacy_stem)))
        for name in sorted(stn):
            a, b = stn[name], leg[name]
            if a.shape != b.shape:
                sys.exit("{}: shape mismatch {} vs {}".format(
                    name, a.shape, b.shape))
            dev = np.abs(a - b)
            scale = max(np.abs(a).max(), np.abs(b).max())
            rel = dev.max() / scale if scale > 0 else dev.max()
            worst = max(worst, rel)
            print("{:24s} shape {} maxRelErr {:.3e}".format(
                name, a.shape, rel))

    print("worst relative deviation: {:.3e} (tolerance {:.1e})".format(
        worst, args.tol))
    if worst > args.tol:
        sys.exit("FAIL: stencil CB meson fields deviate from legacy "
                 "beyond tolerance")
    print("PASS: StagA2AMesonField (stencil, cbPairs) reproduces "
          "StagA2AMesonFieldLegacy (action) on checkerboarded low modes")


if __name__ == "__main__":
    main()
