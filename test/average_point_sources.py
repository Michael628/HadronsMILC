#!/usr/bin/env python3
"""Average point-source meson correlators over all point sources in a directory.

Reads every corr_*.h5 file in test/point-source-spintaste-reference/, averages
the 'corr' attribute of each of the 48 meson groups over all point source
files, and writes the result to all-points.20.h5 in the same directory.
"""

import glob
import os
import re

import h5py
import numpy as np

SRC_DIR = os.path.join(os.path.dirname(__file__), "point-source-spintaste-reference")
OUT_FILE = os.path.join(SRC_DIR, "all-points.20.h5")
N_MESON = 48
T_RE = re.compile(r"t(\d+)")


def source_t(path):
    m = T_RE.search(os.path.basename(path))
    if not m:
        raise ValueError(f"could not parse source time from {path}")
    return int(m.group(1))


def main():
    files = sorted(glob.glob(os.path.join(SRC_DIR, "corr_*.h5")))
    if not files:
        raise RuntimeError(f"no point source files found in {SRC_DIR}")

    sums = [None] * N_MESON
    meta = [None] * N_MESON

    for path in files:
        t = source_t(path)
        with h5py.File(path, "r") as f:
            for i in range(N_MESON):
                g = f[f"meson/meson_{i}"]
                corr = np.roll(g.attrs["corr"], -t)
                if sums[i] is None:
                    sums[i] = corr.copy().astype(corr.dtype)
                    meta[i] = {
                        "scaling": g.attrs["scaling"],
                        "sinkGamma": g.attrs["sinkGamma"],
                        "sourceGamma": g.attrs["sourceGamma"],
                    }
                else:
                    sums[i]["re"] += corr["re"]
                    sums[i]["im"] += corr["im"]

    n = len(files)
    with h5py.File(OUT_FILE, "w") as out:
        for i in range(N_MESON):
            g = out.create_group(f"meson/meson_{i}")
            g.create_group("timeShifts")

            avg = sums[i].copy()
            avg["re"] /= n
            avg["im"] /= n

            g.attrs["corr"] = avg
            g.attrs["scaling"] = meta[i]["scaling"]
            g.attrs["sinkGamma"] = meta[i]["sinkGamma"]
            g.attrs["sourceGamma"] = meta[i]["sourceGamma"]

    print(f"averaged {n} point source files -> {OUT_FILE}")


if __name__ == "__main__":
    main()
