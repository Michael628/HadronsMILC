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
       --with-lime=${TOPDIR}/deps/install${BUILD_EXT} \
       CXX="nvcc" MPICXX="mpicxx" \
	     LDFLAGS="-cudart shared" \
       CXXFLAGS="-ccbin mpicxx -gencode arch=compute_80,code=sm_80 -std=c++17 -cudart shared" 
}

