#!/usr/bin/env python3
"""Validate the SIB HVP A2A batch (RandomWall `_vec`, StagLMAMesonFieldProp
`a2a_batch`, StagA2AMesonField single-key unwrap) from ONE run of
params/sib-hvp-a2a-batch.20.xml (tStep=2 -> nSlices=2 slices t=0,2; nSrc=2
-> 12 columns in order 3*(noise*nSlices + slice) + color).

Five checks, all against probe meson fields MF(eigfull_vec, .) (identity
gamma G1_G1, zero momentum: each column's entries are its projections onto
the |e+o>/|e-o> pair basis, so column-for-column comparisons are exact
statements about the vector sets):

  1. column-order match : MF(h) ~= MF(g) per column (the low-mode part of
     the unprojected solve IS the guess; equality up to the 1e-8 solve
     residual) -- the bit-match gate for the whole column-order contract.
  2. batch vs per-slice : MF(g) equals the default-mode per-slice outputs
     (GammaMapElement -> PropToFermions bridge, noise-major within each
     slice) reordered into batch order -- EXACT (same reconstruction path).
  3. RandomWall `_vec`  : MF(`_vec`) column 3*(n*nSlices+j)+c at timeslice
     t equals MF(eta_fv) column 3*n+c at t == t0+j*tStep and zero at every
     other t -- EXACT (where-masking copies the on-slice values).
  4. unwrap equivalence : MF(right=batch map) bitwise equals MF(right=
     GammaMapElement bare copy of the same map).
  5. guess quality      : from the run log, the CG `Outer iteration 0
     residual` of quark_batchguess is no worse than quark_zeroguess's,
     field by field (compare_guess_quality.py's regression gate).

Usage (from test/, after
../build-scalar/HadronsMILC params/sib-hvp-a2a-batch.20.xml --grid 4.4.4.4 \
    > work/sib-a2a-batch.log):

    python3 compare_sib_hvp_a2a_batch.py [--traj 20] [--tol 1e-4] \
        [--tol-exact 0.0] [--log work/sib-a2a-batch.log]
"""

import argparse
import os
import re
import sys

import h5py
import numpy as np

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATASET = "G1_G1_0_0_0"

NSRC = 2
NSLICES = 2
T0 = 0
TSTEP = 2

MODULE_RE = re.compile(r"module '([^']+)'\)")
RESIDUAL_RE = re.compile(r"Outer iteration 0 residual (\S+) target")


def load(stem, traj):
    path = os.path.join(TEST_DIR, "{}.{}".format(stem, traj),
                        DATASET + ".h5")
    with h5py.File(path, "r") as f:
        d = f[DATASET]["a2aMatrix"][()]
    return (d["re"] + 1j * d["im"]).astype(np.complex64)


def batch_from_perslice(per_slice, nsrc, nslices):
    """Reorder per-slice noise-major blocks [.., n*3+c] into batch
    column order 3*(n*nslices + j) + c."""
    nt, rows, _ = per_slice[0].shape
    out = np.zeros((nt, rows, nsrc * nslices * 3), dtype=per_slice[0].dtype)
    for j in range(nslices):
        src = per_slice[j]
        for n in range(nsrc):
            for c in range(3):
                out[:, :, 3 * (n * nslices + j) + c] = src[:, :, 3 * n + c]
    return out


def check(name, ok, detail):
    print("[{}] {}: {}".format("PASS" if ok else "FAIL", name, detail))
    return ok


def outer0_residuals(lines, module):
    residuals = []
    current = None
    for line in lines:
        m = MODULE_RE.search(line)
        if m:
            current = m.group(1)
            continue
        if current == module:
            r = RESIDUAL_RE.search(line)
            if r:
                residuals.append(float(r.group(1)))
    return residuals


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", default="work/sib-a2a")
    p.add_argument("--traj", type=int, default=20)
    p.add_argument("--tol", type=float, default=1e-4,
                   help="check 1 relative tolerance (solve residual)")
    p.add_argument("--tol-exact", type=float, default=0.0,
                   help="checks 2-4 tolerance (0 = bitwise)")
    p.add_argument("--log", default="work/sib-a2a-batch.log")
    p.add_argument("--slack", type=float, default=1e-6,
                   help="check 5 multiplicative slack (guess-quality gate)")
    args = p.parse_args()

    def stem(name):
        return os.path.join(args.work, name)

    eta = load(stem("eta"), args.traj)
    rwvec = load(stem("rwvec"), args.traj)
    h = load(stem("h"), args.traj)
    g = load(stem("g"), args.traj)
    gps = [load(stem("gps0"), args.traj), load(stem("gps2"), args.traj)]
    gbare = load(stem("gbare"), args.traj)

    ncol = NSRC * NSLICES * 3
    for name, arr in [("h", h), ("g", g), ("rwvec", rwvec),
                      ("gbare", gbare)]:
        if arr.shape[2] != ncol:
            sys.exit("FAIL: {} has {} columns, expected {}".format(
                name, arr.shape[2], ncol))
    for j, arr in enumerate(gps):
        if arr.shape[2] != NSRC * 3:
            sys.exit("FAIL: gps{} has {} columns, expected {}".format(
                j, arr.shape[2], NSRC * 3))

    ok = True

    # 4. unwrap equivalence (bitwise)
    dev = np.abs(g.astype(np.complex128) - gbare.astype(np.complex128)).max()
    ok &= check("4 unwrap equivalence",
                dev <= args.tol_exact, "max |delta| = {}".format(dev))

    # 2. batch vs per-slice (exact)
    g_ps = batch_from_perslice(gps, NSRC, NSLICES)
    dev = np.abs(g.astype(np.complex128) -
                 g_ps.astype(np.complex128)).max()
    scale = np.abs(g).max()
    rel = dev / scale if scale > 0 else dev
    ok &= check("2 batch vs per-slice", rel <= args.tol_exact,
                "max |delta| = {} (rel {:.3e})".format(dev, rel))

    # 1. column-order match: low-mode part of the solve == the guess
    dev = np.abs(h.astype(np.complex128) - g.astype(np.complex128))
    scale = np.abs(g).max()
    rel = dev.max() / scale if scale > 0 else dev.max()
    ok &= check("1 column-order match", rel <= args.tol,
                "max relative deviation {:.3e} (tol {:.1e})".format(
                    rel, args.tol))

    # 3. RandomWall _vec: masked eta per (noise, slice, color)
    worst = 0.0
    exact = True
    for t in range(eta.shape[0]):
        for n in range(NSRC):
            for j in range(NSLICES):
                tj = T0 + j * TSTEP
                for c in range(3):
                    col = 3 * (n * NSLICES + j) + c
                    expect = (eta[t, :, 3 * n + c] if t == tj
                              else np.zeros(eta.shape[1],
                                            dtype=eta.dtype))
                    d = np.abs(
                        rwvec[t, :, col].astype(np.complex128) -
                        expect.astype(np.complex128)).max()
                    worst = max(worst, d)
                    exact &= (d <= args.tol_exact)
    ok &= check("3 RandomWall _vec", exact,
                "worst |delta| = {} (masked slices must vanish)".format(
                    worst))

    # 5. guess quality from the run log
    with open(os.path.join(TEST_DIR, args.log)) as f:
        lines = f.readlines()
    zero = outer0_residuals(lines, "quark_zeroguess")
    batch = outer0_residuals(lines, "quark_batchguess")
    if not zero or len(batch) != len(zero):
        ok &= check("5 guess quality", False,
                    "log parse: zero {} / batch {} residuals".format(
                        len(zero), len(batch)))
    else:
        ratios = [fb / fz if fz > 0 else float("inf")
                  for fz, fb in zip(zero, batch)]
        worst_r = max(ratios)
        no_worse = all(fb <= fz * (1.0 + args.slack)
                       for fz, fb in zip(zero, batch))
        ok &= check("5 guess quality", no_worse,
                    "{} fields, worst batch/zero ratio {:.3f} "
                    "(strictly < 1 expected with production eig "
                    "counts)".format(len(zero), worst_r))

    print("")
    if ok:
        print("PASS: all five SIB HVP A2A batch checks")
    else:
        sys.exit("FAIL: at least one check failed")


if __name__ == "__main__":
    main()
