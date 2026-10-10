#!/usr/bin/env python3
"""Negative-case harness: assert that a HadronsMILC parameter file
ABORTS at setup with the expected HADRONS_ERROR message (the first
negative-test pattern in test/; every compare_*.py is a positive HDF5
diff). Runs the binary, captures stdout+stderr, and requires BOTH a
non-zero exit AND the expected substring. With no arguments, runs all
three known SIB HVP negative cases.

The abort fires in the scheduler's setup dry-run pass (before any
module executes -- no eigensolver, no gauge load), so each case takes
seconds. The two StagLMAMesonFieldProp cases reference the Phase 4
main run's meson-field file in their loaders; that path is never read
at the dry-run stage, but run params/sib-hvp-a2a-batch.20.xml first if
a scheduler without the dry-run pass is ever used.

Usage (from test/):
    python3 check_setup_error.py                      # all three cases
    python3 check_setup_error.py --params params/sib-hvp-neg-multimap.20.xml \\
        --expect "accepts at most one"
"""

import argparse
import os
import subprocess
import sys

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

CASES = {
    "params/sib-hvp-neg-multimap.20.xml": "accepts at most one",
    "params/sib-hvp-neg-twolabels.20.xml": "requires exactly one label",
    "params/sib-hvp-neg-projector.20.xml": "incompatible with",
}


def run_case(binary, params, expect, grid):
    if not os.path.exists(os.path.join(TEST_DIR, params)):
        return False, "params file not found"
    # resolve the binary against TEST_DIR (not the caller's cwd) so the
    # script works from any directory -- the default
    # ../build-scalar/HadronsMILC is relative to test/, matching the
    # params/cwd anchoring above (Step-8 code-review finding)
    binpath = (binary if os.path.isabs(binary)
               else os.path.join(TEST_DIR, binary))
    proc = subprocess.run([binpath, params, "--grid", grid],
                          cwd=TEST_DIR, capture_output=True, text=True)
    out = proc.stdout + proc.stderr
    hit = expect in out
    aborted = proc.returncode != 0
    detail = "exit {} expect '{}' {}".format(
        proc.returncode, expect, "FOUND" if hit else "NOT FOUND")
    if aborted and hit:
        # show the error line for the run log
        for line in out.splitlines():
            if expect in line:
                print("      | " + line.strip())
                break
    return (aborted and hit), detail


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--params", default=None,
                   help="single params file (requires --expect)")
    p.add_argument("--expect", default=None,
                   help="expected HADRONS_ERROR substring")
    p.add_argument("--binary", default="../build-scalar/HadronsMILC")
    p.add_argument("--grid", default="4.4.4.4")
    args = p.parse_args()

    if args.params:
        if not args.expect:
            sys.exit("FAIL: --expect is required with --params")
        cases = [(args.params, args.expect)]
    else:
        cases = sorted(CASES.items())

    ok = True
    for params, expect in cases:
        passed, detail = run_case(args.binary, params, expect, args.grid)
        print("[{}] {}: {}".format("PASS" if passed else "FAIL",
                                   params, detail))
        ok &= passed

    print("")
    if not ok:
        sys.exit("FAIL: at least one negative case did not abort with "
                 "its expected setup error")
    print("PASS: all negative cases aborted at setup with the "
          "expected Argument errors")


if __name__ == "__main__":
    main()
