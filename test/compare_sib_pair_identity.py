#!/usr/bin/env python3
"""Pair-basis identity harness: the SIB derivable blocks are offline
arithmetic from stored files (research 2026-10-04_17-36-40; the pyfm
producer no longer emits lp/np/scalar-nl). From ONE run of
params/sib-pair-identity.20.xml (4^4; t0=0/tStep=2/nSrc=2 -> nSlices=2,
12 batch columns 3*(n*2+j)+c; epack 20 pairs; shift 0.01), five checks:

  1. probe = w.tab (scalar)   : time-SUMMED eigfull probe rows equal
     C * (w o tab[s_col]) per batch column, C a single global constant
     (C ~ norm0, |C| ~ 1, arg ~ 0). Per-timeslice probe rows carry the
     slice-restricted Gram matrix, NOT the identity -- only the sum
     over t diagonalizes it (the operator-level statement).
  2. lp = ll.w.tab (absolute): per timeslice, per family (scalar
     G1_G1 + vec_local GX_GX/GY_GY/GZ_GZ),
        lp[t][:, col(n,j,c)] = (ll[t] * uf) @ (w o tab[t0+j*tStep])[:, 3n+c]
     where uf un-weights ll's ket-side eval folding
     (uf[2a] = eval_a, uf[2a+1] = conj(eval_a); A2AMatrix.hpp:685-687
     folding, LoadMesonField.hpp:388-397 documented inverse). Norm
     ladder cancels exactly at pairScale = sqrt(2) (production default).
  3. np = tab^dag.w.tab (absolute): per timeslice,
        np[t][:, col(n,j,c)] = conj(tab[t]).T @ (w o tab[s_col])[:, 3n+c]
     -- tab-only; per-slice addressing (bra side reads tab[t], the
     weight side reads the SOURCE slice s_col).
  4. nl_s = eval-weighted conj(tab)^T (absolute):
        nl[t][3n'+c', 2a]   = conj(tab[t][2a, 3n'+c']) / eval_a
        nl[t][3n'+c', 2a+1] = conj(tab[t][2a+1, 3n'+c']) / conj(eval_a)
  5. <p|p> closure (ratio-constant): time-summed diag(pp) per column =
     D * sum_a invmag_a * (|tab[s][2a,col]|^2 + |tab[s][2a+1,col]|^2),
     D a single global constant (D ~ norm0^2/2, not predictable
     offline -- constancy IS the identity).

Weights from the BASE-pack evals file (lambda^2, RealD):
  lam_k = sqrt(evals[k]); m = 2*mass (MILC convention -- must match the
  <mass> of the eig ModifyEigenPackMILC module in the params XML,
  default 0.01); invmag_k = 1/(m^2 + lam_k^2);
  w acts on stored row pairs (2k, 2k+1) DIAGONALLY:
  row 2k *= invmag_k*(m - i lam_k), row 2k+1 *= invmag_k*(m + i lam_k).
  The stored rows are the (e+-o)/sqrt(2) rotation of the research
  artifact's (e_cb, e_!cb) doublet, in which invmag*(m sigma0 - i lam
  sigma_x) is diagonal (sigma_x -> sigma_z under that rotation).
  Empirically pinned by validation 2026-10-04_22-45-14: a per-pair 2x2
  fit of the probe against tab gives off-diagonals < 1e-8 and diagonal
  entries invmag*(m -/+ i lam); all five identities then close at
  <= 5.5e-8 with the probe constant C = 1.
All row-pair components (both parity channels) participate in every
check; the worst per-pair deviation is reported per row component.
Expected: max rel dev ~6e-8 (ComplexF storage) at --tol 1e-6.

Usage (from test/, after
../build-scalar/HadronsMILC params/sib-pair-identity.20.xml --grid 4.4.4.4
    > work/sib-pair-identity.log):

    python3 compare_sib_pair_identity.py [--work work/sib-pair-identity] \
        [--traj 20] [--mass 0.01] [--tol 1e-6] [--atol 1e-12]
"""

import argparse
import os
import sys

import h5py
import numpy as np

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATASET = "G1_G1_0_0_0"
VLOC_GAMMAS = ("GX_GX", "GY_GY", "GZ_GZ")

NSRC = 2
NSLICES = 2
T0 = 0
TSTEP = 2


def load(stem, traj, dataset=DATASET):
    path = os.path.join(TEST_DIR, "{}.{}".format(stem, traj),
                        dataset + ".h5")
    with h5py.File(path, "r") as f:
        d = f[dataset]["a2aMatrix"][()]
    # compound {re, im} storage -> complex view before any arithmetic
    return (d["re"] + 1j * d["im"]).astype(np.complex128)


def load_lam(work, traj):
    # saveResult flat file: <work>/evals.<traj>.h5 -- EigenPackExtractEvals
    # calls saveResult(output, "evals", ...), and Module.cpp resultFilename
    # joins output dir + name + "." + traj + ".h5". NOT a per-traj directory
    # like the meson fields, and NOT <work>.<traj>.h5.
    path = os.path.join(TEST_DIR, work, "evals.{}.h5".format(traj))
    with h5py.File(path, "r") as f:
        ev = np.array(f["evals"][()], dtype=np.float64)
    if ev.ndim != 1 or ev.size < 1:
        sys.exit("evals file {}: expected a 1-d RealD dataset, got {}".format(
            path, ev.shape))
    return np.sqrt(ev)


def w_rows(lam, mass):
    """Per-row pair weight, diagonal in the stored row basis: row 2k
    gets invmag*(m - i*lam), row 2k+1 invmag*(m + i*lam). The stored
    rows are the (e+-o)/sqrt(2) rotation of the (e_cb, e_!cb) doublet
    in which the research's invmag*(m*sigma0 - i*lam*sigma_x) acts
    diagonally (sigma_x -> sigma_z under that rotation)."""
    m = 2.0 * mass
    invmag = 1.0 / (m * m + lam * lam)
    w = np.empty(2 * lam.shape[0], dtype=np.complex128)
    w[0::2] = invmag * (m - 1j * lam)
    w[1::2] = invmag * (m + 1j * lam)
    return w


def ket_unweight(lam, mass):
    """Inverse of the writer's ket-side eval folding (LoadMesonField
    form): stored col 2a *= eval_a, col 2a+1 *= conj(eval_a)."""
    m = 2.0 * mass
    ev = m + 1j * lam
    uf = np.empty(2 * lam.shape[0], dtype=np.complex128)
    uf[0::2] = ev
    uf[1::2] = np.conjugate(ev)
    return uf


def col_decode(col):
    """Batch column 3*(n*NSLICES + j) + c -> (n, source slice, 3n+c)."""
    n, rem = divmod(col, 3 * NSLICES)
    j, c = divmod(rem, 3)
    return n, T0 + j * TSTEP, 3 * n + c


def check(name, ok, detail):
    print("[{}] {}: {}".format("PASS" if ok else "FAIL", name, detail))
    return ok


def rel_dev(dev, ref):
    scale = np.abs(ref).max() if ref.size else 0.0
    return float(dev.max() / scale) if scale > 0 else float(dev.max())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", default="work/sib-pair-identity")
    p.add_argument("--traj", type=int, default=20)
    p.add_argument("--mass", type=float, default=0.01,
                   help="the ModifyEigenPackMILC shift of the params XML "
                        "(the m of w is 2*mass, MILC convention)")
    p.add_argument("--tol", type=float, default=1e-6)
    p.add_argument("--atol", type=float, default=1e-12,
                   help="check 4 absolute floor (atol + tol*|nl|)")
    args = p.parse_args()

    def stem(name):
        return os.path.join(args.work, name)

    lam = load_lam(args.work, args.traj)
    n_eig = lam.shape[0]
    wrow = w_rows(lam, args.mass)
    uf = ket_unweight(lam, args.mass)

    tab = load(stem("tab"), args.traj)                    # (nt, 2nE, 3*NSRC)
    ll = load(stem("s_ll"), args.traj)                    # (nt, 2nE, 2nE)
    lp = load(stem("s_lp"), args.traj)                    # (nt, 2nE, 12)
    np_blk = load(stem("s_np"), args.traj)                # (nt, 3*NSRC, 12)
    nl = load(stem("s_nl"), args.traj)                    # (nt, 3*NSRC, 2nE)
    probe = load(stem("probe_g"), args.traj)              # (nt, 2nE, 12)
    pp = load(stem("pp"), args.traj)                      # (nt, 12, 12)
    vl_ll = [load(stem("vl_ll"), args.traj, g + "_0_0_0") for g in VLOC_GAMMAS]
    vl_lp = [load(stem("vl_lp"), args.traj, g + "_0_0_0") for g in VLOC_GAMMAS]

    nt = tab.shape[0]
    ncol = NSRC * NSLICES * 3
    # loud shape guards (never silently truncate -- 3f6faad lesson)
    for name, arr, shape in (
        ("tab", tab, (nt, 2 * n_eig, 3 * NSRC)),
        ("s_ll", ll, (nt, 2 * n_eig, 2 * n_eig)),
        ("s_lp", lp, (nt, 2 * n_eig, ncol)),
        ("s_np", np_blk, (nt, 3 * NSRC, ncol)),
        ("s_nl", nl, (nt, 3 * NSRC, 2 * n_eig)),
        ("probe_g", probe, (nt, 2 * n_eig, ncol)),
        ("pp", pp, (nt, ncol, ncol)),
    ):
        if arr.shape != shape:
            sys.exit("FAIL: {} has shape {}, expected {}".format(
                name, arr.shape, shape))

    # per-source-slice weighted tabs: wt[j] = w o tab[t0+j*tStep]
    # (diagonal row scaling)
    wt = [tab[T0 + j * TSTEP] * wrow[:, np.newaxis] for j in range(NSLICES)]

    ok = True

    # ---- check 2: lp = (ll*uf) . (w o tab[s_col]) per timeslice ----
    def lp_check(name, ll_arr, lp_arr):
        ll_unw = ll_arr * uf[np.newaxis, np.newaxis, :]
        derived = np.empty_like(lp_arr)
        for col in range(ncol):
            _, s, tc = col_decode(col)
            derived[:, :, col] = np.matmul(
                ll_unw, wt[(s - T0) // TSTEP][:, tc])
        dev = np.abs(lp_arr - derived)
        worst = worst_even = worst_odd = 0.0
        worst_pair = -1
        for a in range(n_eig):
            rows = slice(2 * a, 2 * a + 2)
            scale = np.abs(lp_arr[:, rows, :]).max()
            if scale <= 0:
                continue
            rel = float(dev[:, rows, :].max() / scale)
            if rel > worst:
                worst, worst_pair = rel, a
            worst_even = max(worst_even,
                             float(dev[:, 2 * a, :].max() / scale))
            worst_odd = max(worst_odd,
                            float(dev[:, 2 * a + 1, :].max() / scale))
        return check(
            "2 lp = ll.w.tab ({})".format(name),
            worst <= args.tol,
            "worst per-pair rel dev {:.3e} at pair {} "
            "(row 2k {:.3e} / row 2k+1 {:.3e}; tol {:.1e})".format(
                worst, worst_pair, worst_even, worst_odd, args.tol))

    ok &= lp_check("G1_G1", ll, lp)
    for g, a_ll, a_lp in zip(VLOC_GAMMAS, vl_ll, vl_lp):
        ok &= lp_check(g, a_ll, a_lp)

    # ---- check 3: np = conj(tab[t]).T . (w o tab[s_col]) ----
    derived_np = np.empty_like(np_blk)
    for t in range(nt):
        for col in range(ncol):
            _, s, tc = col_decode(col)
            derived_np[t, :, col] = np.matmul(
                np.conjugate(tab[t]).T, wt[(s - T0) // TSTEP][:, tc])
    rel3 = rel_dev(np.abs(np_blk - derived_np), np_blk)
    ok &= check("3 np = tab^dag.w.tab", rel3 <= args.tol,
                "max rel dev {:.3e} (tol {:.1e})".format(rel3, args.tol))

    # ---- check 4: nl_s = eval-weighted conj(tab)^T ----
    ev = 2.0 * args.mass + 1j * lam
    nl_derived = np.empty_like(nl)
    nl_derived[:, :, 0::2] = np.conjugate(
        np.transpose(tab, (0, 2, 1))[:, :, 0::2]) / ev[np.newaxis, np.newaxis, :]
    nl_derived[:, :, 1::2] = np.conjugate(
        np.transpose(tab, (0, 2, 1))[:, :, 1::2]) / np.conjugate(
            ev[np.newaxis, np.newaxis, :])
    dev4 = np.abs(nl - nl_derived)
    within4 = bool(np.all(dev4 <= args.atol + args.tol * np.abs(nl)))
    ok &= check("4 nl_s = eval-weighted conj(tab)^T", within4,
                "max |delta| {:.3e} (atol {:.1e} + tol {:.1e}*|nl|)".format(
                    dev4.max(), args.atol, args.tol))

    # ---- check 1: time-summed probe = C * (w o tab[s_col]) ----
    probe_tsum = probe.sum(axis=0)
    rhs = np.empty_like(probe_tsum)
    for col in range(ncol):
        _, s, tc = col_decode(col)
        rhs[:, col] = wt[(s - T0) // TSTEP][:, tc]
    rnorm = float(np.vdot(rhs.ravel(), rhs.ravel()).real)
    if rnorm <= 0:
        sys.exit("FAIL: probe right-hand side vanishes")
    c_const = np.vdot(rhs.ravel(), probe_tsum.ravel()) / rnorm
    resid = np.abs(probe_tsum - c_const * rhs)
    rel1 = rel_dev(resid, probe_tsum)
    ok &= check("1 probe = C*(w o tab) [time-summed]", rel1 <= args.tol,
                "residual rel dev {:.3e}; C = {:.6f}{:+.6f}j "
                "(expect |C| ~ 1, arg ~ 0)".format(
                    rel1, c_const.real, c_const.imag))

    # ---- check 5: <p|p> closure (ratio-constant) ----
    diag = np.array([pp[t][c, c] for t in range(nt)
                     for c in range(ncol)]).reshape(nt, ncol).sum(axis=0)
    invmag = 1.0 / ((2.0 * args.mass) ** 2 + lam ** 2)
    tnorm = np.empty(ncol)
    for col in range(ncol):
        _, s, tc = col_decode(col)
        x = tab[s][:, tc].reshape(n_eig, 2)
        tnorm[col] = float(np.sum(invmag * (np.abs(x) ** 2).sum(axis=1)))
    if np.any(tnorm <= 0):
        sys.exit("FAIL: <p|p> tab-side norm vanishes for some column")
    ratio = diag / tnorm
    d_mean = ratio.mean()
    spread = float(np.abs(ratio - d_mean).max() / abs(d_mean))
    ok &= check("5 <p|p> closure", spread <= args.tol,
                "ratio spread {:.3e} about D = {:.6f} "
                "(D ~ norm0^2/2; constancy is the identity)".format(
                    spread, d_mean))

    print("")
    if ok:
        print("PASS: all five pair-basis identity checks")
    else:
        sys.exit("FAIL: at least one check failed")


if __name__ == "__main__":
    main()
