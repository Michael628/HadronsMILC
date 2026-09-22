#!/usr/bin/env python3
"""Bit-compare the file-driven LMA propagator producer
(MFermion::StagLMAMesonFieldProp) against the live LowModeProj
(MSolver::StagLMA, projector branch) on the 4^4 test configuration.

The producer reconstructs one PropagatorField per noise per timeslice
per gamma from the loaded meson-field tables (nNoise windows of 3
adjacent color columns; nNoise == 1 collapses to a scalar propagator),
assembling color slot c of noise n from table column
noiseIndex + 3n + c. MUtilities::PropToFermions flattens the output(s)
NOISE-MAJOR into a FermionField vector of 3*nNoise entries, which the
probe contraction (MContraction::StagA2AMesonField, identity spin-taste,
zero momentum) contracts against the full |e+o>/|e-o> eigenvector-pair
basis. The reference side applies the live solver to every noise field
through GaugeProp; its columns col..col+3*nNoise-1 are the per-noise
color components at the same timeslice, so file column 3n + c is
compared against reference column col + 3n + c element by element --
EVERY noise window is checked (conjugation/phase/normalization bugs
have historically surfaced per window, fd87d4f precedent).

Usage (from the test/ directory, after running the schedule
params/lma-mesonfield-file-compare.20.xml with ../HadronsMILC --grid 4.4.4.4):

    python3 compare_lma_mesonfield.py [--traj 20] [--col 0] [--nnoise 1] [--tol 1e-4]
"""

import argparse
import os
import sys

import h5py
import numpy as np

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATASET = "G1_G1_0_0_0"


def load(stem, traj):
    path = os.path.join(TEST_DIR, "{}.{}".format(stem, traj), DATASET + ".h5")
    with h5py.File(path, "r") as f:
        d = f[DATASET]["a2aMatrix"][()]
    # the writer stores an HDF5 compound type (re/im float pair), not a
    # native complex type -- view it as complex64 before any arithmetic
    return (d["re"] + 1j * d["im"]).astype(np.complex64)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ref", default="work/lma-file-compare/ref")
    p.add_argument("--file", default="work/lma-file-compare/file")
    p.add_argument("--traj", type=int, default=20)
    p.add_argument("--col", type=int, default=0,
                   help="base noise column of the color-diluted window "
                        "(noiseIndex); the reference side columns "
                        "col..col+3*nnoise-1 are compared against the "
                        "file side's full noise-major column set")
    p.add_argument("--nnoise", type=int, default=1,
                   help="number of noise windows compared (nNoise); the "
                        "file side holds 3*nnoise noise-major columns "
                        "(3 color columns per window)")
    p.add_argument("--tol", type=float, default=1e-4)
    args = p.parse_args()

    # 3 color columns per noise window (FImpl::Dimension), noise-major
    ncol = 3 * args.nnoise

    ref = load(args.ref, args.traj)[:, :, args.col:args.col + ncol]
    fil = load(args.file, args.traj)
    # fail loudly on an nNoise mismatch instead of silently truncating
    # (the file side is the FULL noise-major column set: 3 per window)
    if fil.shape[2] != ncol:
        sys.exit("file side has {} columns, expected 3*nnoise = {} -- "
                 "did the run use a different nNoise?".format(
                     fil.shape[2], ncol))

    if ref.shape != fil.shape:
        sys.exit("shape mismatch: {} vs {}".format(ref.shape, fil.shape))

    dev = np.abs(fil - ref)
    scale = np.abs(ref).max()
    rel = dev.max() / scale if scale > 0 else dev.max()
    nonzero = dev[dev > 0]
    med = np.median(nonzero) if nonzero.size else 0.0

    print("probe contraction shape : {} ({} noise window(s) x 3 color "
          "columns; reference columns {}..{})".format(
              ref.shape, args.nnoise, args.col, args.col + ncol - 1))
    print("max |delta|            : {:.3e}".format(dev.max()))
    print("median nonzero |delta| : {:.3e}".format(med))
    print("max relative deviation : {:.3e} (tolerance {:.1e})".format(rel, args.tol))

    if rel > args.tol:
        sys.exit("FAIL: relative deviation exceeds tolerance")
    print("PASS: StagLMAMesonFieldProp reproduces StagLMA (projector branch)")


if __name__ == "__main__":
    main()
