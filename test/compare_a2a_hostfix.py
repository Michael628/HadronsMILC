#!/usr/bin/env python3
"""Bitwise-compare A2A meson-field HDF5 output trees from two builds of
the same schedule -- the validation gate for the A2A host-loops round
(per-block instrumentation, streaming per-rank conversion, dead-zero
deletion, reconstructLegacy device fusion), whose correctness bar is
bit-identical a2aMatrix datasets on single-rank 4^4 runs.

Both trees come from running the SAME generated params-XML pair (see
work/gen-a2a-hostfix-xml.py) with the golden and fix builds; identical
seeds and inputs make every file value-comparable. The writer schema is
one file per (e,s): top-level group <basename> holding metadata and the
a2aMatrix dataset (nt, ni, nj) stored as a re/im float32 compound (view
as complex64 before arithmetic -- compare_cb_mesonfield.py precedent).
Value-level comparison, never file hashes (chunking may differ); files
without an a2aMatrix dataset (eigenvector packs written by the same
schedule) are skipped with a note.

Usage (from test/, after both runs):

    python3 compare_a2a_hostfix.py [--tol-exact 0.0] \
        work/a2a-hostfix-gt work/a2a-hostfix-devfix

Exit status: 0 iff every a2aMatrix file matches within --tol-exact and
the file sets are equal.
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np


def load(path):
    """Return the a2aMatrix dataset of one file as complex64, or None."""
    with h5py.File(path, "r") as f:
        for key in f.keys():
            if isinstance(f[key], h5py.Group) and "a2aMatrix" in f[key]:
                d = f[key]["a2aMatrix"][()]
                # HDF5 compound (re/im float pair), not a native complex
                # type -- view as complex64 before any arithmetic
                return (d["re"] + 1j * d["im"]).astype(np.complex64)
    return None


def collect(tree):
    """Map relative .h5 path -> complex64 a2aMatrix for one tree."""
    files = sorted(glob.glob(os.path.join(tree, "**", "*.h5"),
                             recursive=True))
    if not files:
        sys.exit("no HDF5 files under {}".format(tree))
    mats, skipped = {}, []
    for path in files:
        rel = os.path.relpath(path, tree)
        m = load(path)
        if m is None:
            skipped.append(rel)
        else:
            mats[rel] = m
    for rel in skipped:
        print("skip {} (no a2aMatrix dataset)".format(rel))
    return mats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("golden", help="golden-build output tree root")
    p.add_argument("fix", help="fix-build output tree root")
    p.add_argument("--tol-exact", type=float, default=0.0,
                   help="max |delta| per dataset (default 0 = bitwise)")
    args = p.parse_args()

    a = collect(args.golden)
    b = collect(args.fix)

    ok = True
    # File-set equality first: a mis-built per-rank I/O task list shows
    # up as missing/extra files before any value comparison bites.
    if set(a) != set(b):
        print("[FAIL] file set mismatch: only-golden {} only-fix {}".format(
            sorted(set(a) - set(b)), sorted(set(b) - set(a))))
        ok = False

    worst, worst_rel, n = 0.0, 0.0, 0
    for rel in sorted(set(a) & set(b)):
        x, y = a[rel], b[rel]
        if x.shape != y.shape:
            print("[FAIL] {}: shape mismatch {} vs {}".format(
                rel, x.shape, y.shape))
            ok = False
            continue
        dev = np.abs(x.astype(np.complex128) -
                     y.astype(np.complex128)).max()
        scale = max(np.abs(x).max(), np.abs(y).max())
        rel_dev = dev / scale if scale > 0 else dev
        worst, worst_rel = max(worst, dev), max(worst_rel, rel_dev)
        n += 1
        good = dev <= args.tol_exact
        ok &= good
        print("[{}] {:44s} shape {} max|delta| {:.3e} (rel {:.3e})".format(
            "PASS" if good else "FAIL", rel, x.shape, dev, rel_dev))

    print("")
    print("{} a2aMatrix files compared, worst |delta| {:.3e} "
          "(rel {:.3e}), tolerance {:.1e}".format(
              n, worst, worst_rel, args.tol_exact))
    if not ok:
        sys.exit("FAIL: A2A host-fix output trees deviate")
    print("PASS: fix build reproduces golden build " +
          ("bitwise" if args.tol_exact == 0.0 else "within tolerance"))


if __name__ == "__main__":
    main()
