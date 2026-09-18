#!/usr/bin/env python3
"""Bit-compare the file-driven LMA solver (MSolver::StagLMAMesonField) against
the live LowModeProj (MSolver::StagLMA, projector branch) on the 4^4 test
configuration.

Both solvers are applied to the same time-diluted noise source; their output
fields are contracted against the full |e+o>/|e-o> eigenvector-pair basis
(MContraction::StagA2AMesonField, identity spin-taste, zero momentum) and the
two a2aMatrix HDF5 datasets are compared element by element on the bound
noise column. The file solver is bound to (noiseIndex, timeslice) and returns
the same field regardless of which noise vector it is applied to, so only
column noiseIndex of the probe contraction is comparable; all other columns
of the file-side matrix repeat that field by design. Precedent: fd87d4f
(conjugation issues found by reference comparison).

Usage (from the test/ directory, after running the schedule
params/lma-mesonfield-file-compare.20.xml with ../HadronsMILC --grid 4.4.4.4):

    python3 compare_lma_mesonfield.py [--traj 20] [--col 0] [--tol 1e-4]
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
    p.add_argument("--ref", default="lma-file-compare/ref")
    p.add_argument("--file", default="lma-file-compare/file")
    p.add_argument("--traj", type=int, default=20)
    p.add_argument("--col", type=int, default=0,
                   help="noise column the file solver is bound to "
                        "(noiseIndex); only this column is comparable")
    p.add_argument("--tol", type=float, default=1e-4)
    args = p.parse_args()

    ref = load(args.ref, args.traj)[:, :, [args.col]]
    fil = load(args.file, args.traj)[:, :, [args.col]]

    if ref.shape != fil.shape:
        sys.exit("shape mismatch: {} vs {}".format(ref.shape, fil.shape))

    dev = np.abs(fil - ref)
    scale = np.abs(ref).max()
    rel = dev.max() / scale if scale > 0 else dev.max()
    nonzero = dev[dev > 0]
    med = np.median(nonzero) if nonzero.size else 0.0

    print("probe contraction shape: {} (column {})".format(ref.shape, args.col))
    print("max |delta|            : {:.3e}".format(dev.max()))
    print("median nonzero |delta| : {:.3e}".format(med))
    print("max relative deviation : {:.3e} (tolerance {:.1e})".format(rel, args.tol))

    if rel > args.tol:
        sys.exit("FAIL: relative deviation exceeds tolerance")
    print("PASS: StagLMAMesonField reproduces StagLMA (projector branch)")


if __name__ == "__main__":
    main()
