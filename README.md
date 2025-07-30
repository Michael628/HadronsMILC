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

## Building

### Prerequisites

This project requires several dependencies and the Grid/Hadrons frameworks to be built and installed. A `build.sh` script has been provided to simplify the build process. The syntax for building each of the four build-steps below is `./build.sh <build-type> <step>`, where `<build-type>` in the provided script can be`scalar` or `mpi`. These are example scripts intended for simple local builds of Grid, and will need to be extended for more sophisticated machines. See the systems directory in Grid's github: <https://github.com/paboyle/Grid/tree/develop/systems> for examples of more complex builds. The `<step>` can be `deps`, `grid`, `hadrons`, or `app`, as described below.

### 0. Build Dependencies

First, install required libraries (GMP, MPFR, HDF5, and LIME):

```bash
# Use the provided dependency script
./build.sh <build-type> deps
```

Or manually install the dependencies using your system package manager.

### 1. Build Grid

```bash
# Use the provided build script for MILC-specific Grid branch
./build.sh <build-type> grid
```

### 2. Build Hadrons

```bash
# Use the provided build script for MILC-specific Hadrons branch
./build.sh <build-type> hadrons
```

### 3. Build HadronsMILC

```bash
# Use the provided build script
./build.sh <build-type app
```

## Usage

Run with XML parameter files similar to standard Hadrons applications:

```bash
./HadronsMILC parameters.xml --grid 4.4.4.4
```

The application uses the Hadrons framework's module system to construct and execute computational graphs for lattice QCD calculations.
