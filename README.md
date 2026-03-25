# HadronsMILC

A lattice QCD physics application built on the Grid/Hadrons framework for computing hadron correlation functions using staggered fermions.

## Overview

HadronsMILC extends the Hadrons framework with specialized modules for MILC-style staggered fermion calculations. It provides tools for computing meson correlation functions, all-to-all (A2A) vector contractions, and other hadron physics observables on lattice QCD gauge configurations.

## Key Features

- **Staggered Fermion Actions**: HISQ action with MILC mass convention built from smeared fields (smearing not yet implemented)
- **A2A Contractions**: Efficient all-to-all vector computations for hadron correlation functions
- **Meson Correlators**: Specialized modules for computing meson correlation functions with staggered spin-taste structure
- **Solvers**: Mixed-precision conjugate gradient solvers and eigenvalue computations
- **Noise Sources**: Support for various diluted noise source types

## Physics Modules

- **MAction**: Improved staggered fermion actions (3D and 5D)
- **MContraction**: Meson correlation function calculations and A2A field contractions
- **MSolver**: Linear solvers, A2A vector generation, and eigenvalue computations
- **MSource**: Various source types including random walls and sequential sources
- **MFermion**: Fermion field operations and gauge propagators
- **MGauge**: Gauge field utilities and stochastic electromagnetic functions

## Dependencies

HadronsMILC requires the following to be built and installed:

- GMP, MPFR, HDF5, LIME
- **Grid** (`feature/LMI-develop` branch): <https://github.com/milc-qcd/Grid/tree/feature/LMI-develop>
- **Hadrons** (`feature/LMI-develop` branch): <https://github.com/milc-qcd/Hadrons/tree/feature/LMI-develop>

Build orchestration for all dependencies is managed by pyfm: <https://github.com/Michael628/pyfm>

## Building HadronsMILC

Once Grid and Hadrons are installed and `hadrons-config` is in your PATH:

```bash
./configure.sh && make
```

## Usage

```bash
./HadronsMILC parameters.xml [Grid runtime arguments...]
```

The first argument is an XML parameter file defining the module graph. All subsequent arguments are passed directly to Grid (e.g. `--grid 24.24.24.32`, `--mpi 1.1.1.4`, `--shm 1024`).
