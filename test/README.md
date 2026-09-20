# test/

Validation suite for the HadronsMILC binary. These tests answer one question:
**given this XML, does the binary compute the right numbers?** They need the
compiled app, gauge configs, and minutes to hours of runtime.

They are deliberately *not* pyfm's tests. pyfm asks a different question --
does its generator emit the right XML? -- and answers it in seconds with no
binary, via golden files in `pyfm/test/tasks/in/`. Both tests guard the same
XML module-schema contract from opposite sides; neither is redundant.

## The rules

1. **Track a file iff no run writes to its path.** Anything the binary or a
   script produces goes under `work/`, which is gitignored.
2. **A param file lives here iff it is hand-authored.** Anything pyfm
   generates belongs next to the generator that emits it.

Rule 2 is why `in/`, `schedules/`, `params.yaml` and `hadrons.sh` are gone --
they carried pyfm's runId template `LMI-RW-series-{series}-{eigs}-eigs-{noise}-noise`
and were pyfm output living in the wrong repo.

## Layout

| Dir | Tracked | Contents |
|---|---|---|
| `params/` | yes | Hand-authored param files. Paths are relative to `test/`. |
| `lat/` | yes | Gauge config fixtures. |
| `ref/` | yes | Frozen references. **Never a write target.** |
| `work/` | **no** | All run output. |
| `*.py` | yes | Comparison scripts. |

## Running

From this directory, with `HadronsMILC` built:

```bash
../HadronsMILC params/<name>.xml --grid 4.4.4.4
python3 compare_<name>.py
```

Output lands in `work/`. Comparison scripts read `ref/` and `work/`.

| Script | Compares |
|---|---|
| `compare_lma_mesonfield.py` | `StagLMAMesonField` vs `StagLMA` projector branch |
| `compare_cb_mesonfield.py` | stencil CB meson field vs `StagA2AMesonFieldLegacy` |
| `compare_point_source_a2a.py` | A2A LL correlators vs the point-source reference |
| `average_point_sources.py` | averages 256 point sources -> `all-points.20.h5` |

## Promoting a reference

`ref/` is immutable to runs. To update one, regenerate into `work/`, verify,
then copy it across deliberately:

```bash
python3 average_point_sources.py          # writes work/point-source-spintaste/
cp work/point-source-spintaste/all-points.20.h5 ref/point-source-spintaste/
```

Every entry in `ref/` needs a row below saying how it was made. A reference
you cannot regenerate is one you cannot trust when it starts failing.

| Reference | Produced by |
|---|---|
| `ref/point-source-spintaste/all-points.20.h5` | `params/point-source-spintaste-all-sites.xml`, then `average_point_sources.py` |
| `ref/random-wall-spintaste/corr_32_rw.20.h5` | `params/random-wall-spintaste-32.xml` |
| `ref/random-wall-fullvolume/corr_32_rw.20.h5` | `params/random-wall-fullvolume-32.xml` |

## Fixture provenance

| File | Status |
|---|---|
| `lat/lat.sample.l4444.ildg.20` | **live** -- read by most params |
| `lat/fatlinks.l4444.ildg.20` | **live** -- HISQ fat links for the above |
| `lat/longlinks.l4444.ildg.20` | **live** -- HISQ long links for the above |
| `lat/{l,fat,lng}{4444,6666,8888}_free.ildg` | unreferenced; free-field, regenerable |
| `lat/test-{fat,long}links.l4444.ildg.20` | unreferenced |

No param file here loads an eigenpack from disk -- the three that need
eigenvectors generate them with `MSolver::StagFermionIRL`. The former
`eigen/` directory served only the pyfm-generated params and moved out with
them.

## Known gaps

- `ref/random-wall-spintaste/` and `ref/random-wall-fullvolume/` have no
  script that reads them. The naming implies a comparison that was run by
  hand; `work/random-wall-spintaste-regression/` differs from the reference
  and nothing reports it.
- `params/milc-ref-input-point-source-spintaste-src0000.in` is a MILC input,
  not a HadronsMILC one. Its upstream workspace is `../../test/milc/`.
