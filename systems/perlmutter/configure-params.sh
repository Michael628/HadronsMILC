#!/bin/bash
# configure-params.sh - Machine-specific configuration for builds
#
# This script provides functions that set environment variables for configure calls.
# It expects the following BUILD_* variables to be set by build.sh before sourcing:
#   - BUILD_EXT
#   - BUILD_DEBUG
#   - BUILD_MPI_REDUCTION

function grid_configure() {
  local INSTALLDIR=$1
  local TOPDIR=$2
  ${TOPDIR}/Grid/configure \
   --prefix=${INSTALLDIR} \
   --enable-comms=mpi-auto       \
   --enable-simd=GPU \
   --enable-shm=nvlink \
   --enable-gen-simd-width=64 \
   --enable-accelerator=cuda \
   --disable-fermion-reps \
   --disable-unified \
   --disable-gparity \
   --with-mpfr=${TOPDIR}/deps/install${BUILD_EXT} \
   CXX=nvcc \
   LDFLAGS='-cudart shared' \
   CXXFLAGS='-ccbin CC -gencode arch=compute_80,code=sm_80 -I${CUBLAS_PATH}/include -std=c++17 -cudart shared -DEIGEN_DONT_VECTORIZE'  \
   LIBS='-lcublas -lhdf5_cpp -L${CUBLAS_PATH}/lib'
}

function hadrons_configure() {
  local INSTALLDIR=$1
  local TOPDIR=$2

  # Configure arguments for Hadrons
  ${TOPDIR}/Hadrons/configure \
    --prefix=${INSTALLDIR} \
    --with-grid=${TOPDIR}/Grid/install${BUILD_EXT}
}

function app_configure() {
  local INSTALLDIR=$1
  local TOPDIR=$2

  # Configure arguments for App
  ${TOPDIR}/HadronsMILC/configure \
  --prefix=${INSTALLDIR} \
  --with-grid=${TOPDIR}/Grid/install${BUILD_EXT} \
  --with-hadrons=${TOPDIR}/Hadrons/install${BUILD_EXT}

}

function dependency_configure() {
  local dep_name=$1
  local INSTALLDIR=$2

  # Dependency-specific additions
  case ${dep_name} in
    mpfr)
      DEP_CONFIGURE_ARGS="--with-gmp=${INSTALLDIR}"
    ;;
    hdf5)
      DEP_CONFIGURE_ARGS="--enable-cxx"
    ;;
  esac
  ${TOPDIR}/deps/${dep_name}/configure \
    --prefix=${INSTALLDIR} \
    ${DEP_CONFIGURE_ARGS} \
    CXXFLAGS=-O3 \
    CFLAGS=-O3 \
    CC=cc CXX=CC

}
