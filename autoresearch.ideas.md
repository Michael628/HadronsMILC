# Autoresearch ideas — HadronsMILC

## Landed
- `422c7aa` — fuse SUM/DIFF axpy pairs (2 kernels/eigenpair → 1): 18.91 → 9.13 ms.
- `012f1e3` — fuse across color columns (6 accumulators, 3·nEigs → nEigs kernels/output): 8.44 → 7.37 ms. Accumulation is no longer the dominant per-output cost on 4^4.

## Verified-but-discarded (run #9)
- **rbFermNeg zero-init before Meooe**: Meooe DOES fully overwrite its output (ImprovedStaggeredFermion::Meooe → DhopOE/EO → DhopImproved opens out AcceleratorWrite, DhopSiteGeneric coalescedWrites every site; first stencil leg assigns, never reads out) — proven by source AND by bit-identity with stale buffer content. But removing it moved the metric −0.2% (noise floor ~1 ms): discarded. Reapply the one-line removal only if a production-lattice profile shows the per-column tail dominating.

## Open ideas (module ~0.9 ms/output on 4^4; per-column tail + zero-inits dominate)
- **First-iteration write instead of 6 accumulator Zero()s**: the k-loop's first eigenpair (k = eigStart+nEigs−1) could write `coef·ev` unconditionally instead of `coef·ev + read(acc)`, deleting all six zero-init kernels per output. Caveat: `x + 0.0` vs `x` differs on signed zeros — numpy-value gates pass but strict bit-identity would not; decide if the signed-zero caveat is acceptable documentation-wise before trying.
- **outputName()/gammaList() parse hoisting** (re-parses per name; ~2 calls/output in execute()): sub-noise on 4^4, grows with output count; pure hygiene — poor autoresearch candidate, fine as drive-by cleanup.
- **Bigger benchmark needed**: at 4^4 every full-field op is launch-overhead dominated (~8 µs for 6 KB) and box noise (~±1 ms on the module timer) swamps everything below ~10% module-level change. For further hot-loop work, consider a larger-volume scratch schedule (e.g. 8^4 or more evecs via nEigs) so memory-bound effects and kernel-count changes resolve above noise.
- **Grid view-family exclusivity (lesson for custom kernels)**: a buffer's Accelerator-mode views must ALL close before any CPU-mode view (setCheckerboard, reductions, transfers) opens on it — `MemoryManagerCache.cc` asserts. Scope hoisted views tightly.
- **autoView macro pitfall**: `autoView(n, *ptr[i], mode)` expands to `*ptr[i].View(mode)` — `.` binds tighter than `*`; bind a named reference first.
- **Harness pitfalls**: Hadrons output directories (`file.20/`, `ref.20/`) collide across schedules — snapshot per-schedule subdirs; the two schedules' `ref.20` datasets legitimately differ. Box noise is heavy — min-of-3+ per arm and confirm marginal wins with a second full run.
