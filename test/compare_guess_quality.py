#!/usr/bin/env python3
"""Assert that StagLMAMesonFieldProp's file-driven guess actually
deflates MSolver::StagMixedPrecisionCG, using
params/lma-mesonfield-guess-quality.20.xml's three StagGaugePropLegacy
solves (quark_zeroguess / quark_liveguess / quark_fileguess) over the
same noise_vec and the same solver.

Grid's MixedPrecisionConjugateGradient prints one line per outer
iteration per field:

    MixedPrecisionConjugateGradient: Outer iteration 0 residual <r> target <t>

computed AFTER the guess is seeded into the solution vector (see
Grid/algorithms/iterative/ConjugateGradientMixedPrec.h:118-120) -- a
direct, unambiguous measurement of guess quality. This script isolates
each module's "Outer iteration 0" residuals via Hadrons' per-module log
banner ("Measurement step N/M (module '<name>')") and asserts the
file-driven guess's residuals are no worse than the zero-guess
baseline's, field by field -- this is exactly the regression the
missing-eigenvalue-weighting bug caused (file-guess residual grew to
~12x the zero-guess baseline at high eigenvector count).

Usage (from the test/ directory, after running
params/lma-mesonfield-guess-quality.20.xml with
../HadronsMILC --grid 4.4.4.4 > work/guess-quality.log):

    python3 compare_guess_quality.py [--log work/guess-quality.log]
"""

import argparse
import re
import sys

MODULE_RE = re.compile(r"module '([^']+)'\)")
RESIDUAL_RE = re.compile(r"Outer iteration 0 residual (\S+) target")


def outer0_residuals(lines, module):
    """Every 'Outer iteration 0' residual logged during `module`'s
    execute() step, in solve order (one per field)."""
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
    p.add_argument("--log", default="work/guess-quality.log")
    p.add_argument("--zero-module", default="quark_zeroguess")
    p.add_argument("--file-module", default="quark_fileguess")
    p.add_argument("--live-module", default="quark_liveguess")
    p.add_argument("--slack", type=float, default=1e-6,
                   help="multiplicative slack on the zero-guess baseline "
                        "to absorb floating-point noise (default 1e-6)")
    args = p.parse_args()

    with open(args.log) as f:
        lines = f.readlines()

    zero = outer0_residuals(lines, args.zero_module)
    live = outer0_residuals(lines, args.live_module)
    file_ = outer0_residuals(lines, args.file_module)

    if not zero:
        sys.exit("FAIL: no 'Outer iteration 0' residuals found for module "
                 "'{}' -- did the run actually execute "
                 "StagMixedPrecisionCG?".format(args.zero_module))
    if len(file_) != len(zero):
        sys.exit("FAIL: field count mismatch -- zero-guess module '{}' "
                 "logged {} residuals, file-guess module '{}' logged {}"
                 .format(args.zero_module, len(zero), args.file_module,
                         len(file_)))

    print("fields compared         : {}".format(len(zero)))
    for i, (rz, rf) in enumerate(zip(zero, file_)):
        ratio = rf / rz if rz > 0 else float("inf")
        print("  field {}: zero-guess residual {:.6e}, file-guess residual "
              "{:.6e} (ratio {:.3f})".format(i, rz, rf, ratio))

    worst = max((rf / rz if rz > 0 else float("inf"))
               for rz, rf in zip(zero, file_))
    print("worst file/zero residual ratio : {:.3f}".format(worst))
    if live:
        print("(live-guess residuals logged for '{}': {})".format(
            args.live_module, ["{:.6e}".format(r) for r in live]))

    if any(rf > rz * (1.0 + args.slack) for rz, rf in zip(zero, file_)):
        sys.exit("FAIL: file-driven guess produced a WORSE post-guess "
                 "residual than no guess at all for at least one field "
                 "-- this is the missing-eigenvalue-weighting regression")
    print("PASS: StagLMAMesonFieldProp's file-driven guess deflates "
          "StagMixedPrecisionCG (no field regressed past the zero-guess "
          "baseline)")


if __name__ == "__main__":
    main()
