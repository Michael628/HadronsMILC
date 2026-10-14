# Changelog

HadronsMILC shares MAJOR.MINOR with pyfm: `pyfm 0.X.*` works with
`HadronsMILC 0.X.*`. An **Interface** change (a module or parameter added,
renamed or removed, or an output format change) bumps MINOR in both repos.
Repo-internal fixes bump PATCH independently. Dependency pins for each release
are in `DEPENDENCIES`.

## [Unreleased]

### Interface
- `MContraction::StagA2AMesonField` supports general momentum projection:
  any number of non-zero momenta; one HDF5 file per (gamma, momentum) pair
  with the per-file momentum recorded in `md.momentum` (per-file layout and
  the `version = "1"` schema unchanged; zero-momentum outputs unchanged).
- `MContraction::StagA2AASlashMesonField` migrates to the Grid stencil
  worker (J.A local ctor; zero-momentum outputs value-compatible with the
  previous implementation on identical stochastic seeds) and to the
  pair-source CB route: the `action` parameter is REMOVED, replaced by
  `cbPairsLeft`/`cbPairsRight` (two `EigenPackCBPairs` instances, matching
  `StagA2AMesonField`); the module is momentum-capable.
- New `MContraction::StagA2AASlashMesonFieldLegacy`: the frozen
  `A2AWorkerLocal` implementation with the former `action`/swapFn CB route
  (zero-momentum only), retained as a consistency probe for the migration.
- New `MContraction::StagA2AASlashOneLinkMesonField`: per-component
  one-link `J_kA_k` meson fields from raw stochastic photon components
  dressing the gauge links (new required `gauge` parameter; no gammas
  module needed); per-file `amu_index` = the direction k, `amu_seed`
  sample-qualified `<EmSeedString>_<j>`.
- Grid dependency pin bumped to `5325da90c42a61ed76ec407f8156dc662dfe9c8d`
  (general-momentum + J.A local/one-link stencil A2A API; install refresh
  required before building).

## [0.3.0] - 2026-10-10

Requires the GridMilc library and the patched Hadrons pinned in `DEPENDENCIES`.

### Interface
- New `MUtilities` modules: `EigenPackFullPairs`, `EigenPackCBPairs`,
  `EigenPackCheck`, `GammaMapElement` and `PropToFermions`.
- New `MIO` module `LoadMesonField`; new `MFermion` module
  `StagLMAMesonFieldProp`, an eager LMA meson-field propagator producer with
  batched eigenpasses and fused reconstruction (multi-gamma tables, per-noise
  reconstruction, A2A batch output mode, weighted-guess projector parameter).
- The pre-gamma-map implementations are kept as `MFermion::StagGaugePropLegacy`
  and `MContraction::StagMesonLegacy`; `StagGaugeProp` and `StagMeson` are now
  gamma-map-driven.
- `SpinTasteParams` is removed: gamma dispatch in `GaugeProp`, `Meson` and
  `LMAMesonFieldProp` is driven by `SpinTaste`/`GammaMapElement` maps and
  `TGammaMap` outputs.
- `StagA2AMesonField` is stencil-only; `StagA2AMesonFieldOld` is removed.
- `StagRandomWall` supports full-volume external noise.
- `StagA2AVectors` gains SIB comparison parameters and requests a plain solver
  when `lowModes` is given.

### Behavior
- cshift general spin-taste A2A implementation and a stencil-based A2A meson
  field option (performance).
- Reduced memory usage in `StagFermionIRL` setup().
- Meson-field checkerboard guard for low modes.
- Fixed a use-after-free on gauge links.

### Fixed
- `MSolver::StagCGMILC`: the `<name>_subtract` solver returned the guess minus
  the final residual instead of the solution minus the guess. V vectors made
  with `StagA2AVectors` using this solver and guess-vector `lowModes` since the
  module was added are wrong and must be regenerated.

## [0.2.1] - 2026-10-02

Same dependency pins as 0.2.0.

### Added
- `configure` compares the Grid (`Grid/Version.h` `GITHASH`) and Hadrons
  (`hadrons-config --git`) commits against `DEPENDENCIES`. It warns on a
  mismatch, and `--enable-strict-pins` makes a mismatch an error.
- `HadronsMILC --version` and a startup banner report the HadronsMILC version
  and git describe, the Grid and Hadrons commits, the GridMilc version and the
  pin status.
- Optional `<grid><provenance>` element (`pyfmVersion`, `pyfmSha`,
  `hadronsMilcCompat`, `generated`). It is logged when present, and a
  `hadronsMilcCompat` MAJOR.MINOR mismatch produces a warning.

### Changed
- `AC_INIT` version set to `0.2.1`. It said `0.1` through the `v0.2.0` tag.

## [0.2.0] - 2026-10-02

Requires the GridMilc library and the patched Hadrons pinned in `DEPENDENCIES`.

### Interface
- New module `MGauge::HISQSmear` (plus `HISQSmearF` in double-precision
  builds): parameters `gauge` and optional `boundary` (default `1,1,1,-1`);
  outputs `<name>_fat` and `<name>_long`.
- New module `MIO::LoadMilc`: parameters `file` and `exitOnChecksumMismatch`
  (default `false`).
- New module `MIO::SaveIldg`: parameters `gauge`, `fileStem` and optional
  `ensembleLabel`.
- `MAction::ImprovedStaggeredMILC` parameters reduced to `gaugefat`,
  `gaugelong` and `mass`. `c1`, `c2`, `tad`, `boundary`, `string` and `twist`
  were removed. They were already ignored, so there is no physics change.
- Parameter files are parsed by Hadrons' `Application::parseParameterFile`.
  This adds the optional per-module `<subgrid>` tag and
  `<parameters><database><restoreModules>`.

### Behavior
- `MSolver::StagMixedPrecisionCG`: the inner tolerance is now
  `max(residual, 1e-7)`.
- A2A vectors: an empty `lowModes` explicitly means "no low-mode guess".
- `A2AMatrix`: an empty `evals` no longer dereferences `&evals[0]`.

### Build
- `StagGamma` and the A2A worker code moved to GridMilc. The build links
  `-lGridMilc`. XML gamma names are unchanged.

## [0.1.0] - 2026-10-02

Baseline release, the first tagged version. Builds against an unpatched
Hadrons and doesn't need GridMilc.
