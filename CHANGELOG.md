# Changelog

HadronsMILC shares MAJOR.MINOR with pyfm: `pyfm 0.X.*` works with
`HadronsMILC 0.X.*`. An **Interface** change (a module or parameter added,
renamed or removed, or an output format change) bumps MINOR in both repos.
Repo-internal fixes bump PATCH independently. Dependency pins for each release
are in `DEPENDENCIES`.

## [Unreleased]

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
