#!/usr/bin/env python3
"""Compare the averaged point-source meson correlators against the A2A LL
correlators (local + one-link) for the same set of 48 spin-taste operators.

The point-source reference (test/average_point_sources.py output) gives, for
each of the 48 meson operators, a correlator C(dt) already averaged over all
256 point sources and shifted so that index 0 corresponds to the source
timeslice.

The A2A meson field files store C(dt) already averaged over the source time,
one value per separation dt and gamma structure, so they can be compared
directly against the point-source reference.
"""

import os

import h5py
import numpy as np
import pandas as pd

TEST_DIR = os.path.dirname(__file__)
PS_FILE = os.path.join(TEST_DIR, "point-source-spintaste-reference", "all-points.20.h5")
A2A_LOCAL_FILE = os.path.join(
    TEST_DIR,
    "e384n1dt1",
    "correlators",
    "m01",
    "all_local",
    "a2aLL",
    "corr_all_local_a2aLL_m01.20.h5",
)
A2A_ONELINK_FILE = os.path.join(
    TEST_DIR,
    "e384n1dt1",
    "correlators",
    "m01",
    "all_onelink",
    "a2aLL",
    "corr_all_onelink_a2aLL_m01.20.h5",
)
N_MESON = 48
A2A_SCALE = 4**3


def load_point_source_reference(path):
    """Return {gamma: complex ndarray of shape (Lt,)}."""
    result = {}
    with h5py.File(path, "r") as f:
        for i in range(N_MESON):
            g = f[f"meson/meson_{i}"]
            gamma = g.attrs["sinkGamma"][0].decode()
            corr = g.attrs["corr"]
            result[gamma] = corr["re"] + 1j * corr["im"]
    return result


def load_a2a_dt_averaged(path):
    """Return {gamma: complex ndarray of shape (Lt,)} indexed by dt."""
    df = pd.read_hdf(path, key="corr")
    lt = df.index.get_level_values("dt").max() + 1
    result = {}
    for gamma, sub in df.groupby(level="gamma"):
        dt_avg = np.zeros(lt, dtype=complex)
        for (dt, _perm, _gamma), val in sub["corr"].items():
            dt_avg[dt] = val
        result[gamma] = dt_avg / A2A_SCALE
    return result


def main():
    ps = load_point_source_reference(PS_FILE)
    a2a = {}
    a2a.update(load_a2a_dt_averaged(A2A_LOCAL_FILE))
    a2a.update(load_a2a_dt_averaged(A2A_ONELINK_FILE))

    missing = set(ps) - set(a2a)
    if missing:
        print(f"warning: no A2A data found for gammas: {sorted(missing)}")

    header = f"{'gamma':10s} {'dt':>2s} {'point-source':>24s} {'a2a LL':>24s} {'rel diff':>10s} {'sign flip':>9s}"
    print(header)
    print("-" * len(header))

    rel_diffs = []
    for gamma in sorted(ps):
        if gamma not in a2a:
            continue
        ps_corr = ps[gamma]
        a2a_corr = a2a[gamma]
        for dt in range(len(ps_corr)):
            p = ps_corr[dt]
            q = a2a_corr[dt]
            rel = (abs(p) - abs(q)) / abs(p) if abs(p) > 0 else float("nan")
            rel_diffs.append(rel)
            flipped = (p.real * q.real) < 0
            print(f"{gamma:10s} {dt:2d} {p.real!s:>23s} {q.real!s:>24s} {rel:10.4g} {str(flipped):>9s}")

    rel_diffs = np.array(rel_diffs)
    print()
    print(f"max relative diff:  {np.nanmax(rel_diffs):.4g}")
    print(f"mean relative diff: {np.nanmean(rel_diffs):.4g}")


if __name__ == "__main__":
    main()
