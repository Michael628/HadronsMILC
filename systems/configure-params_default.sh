#!/bin/bash
# configure-params.sh - Machine-specific configuration for builds
#
# This script provides functions that set environment variables for configure calls.
# It expects the following BUILD_* variables to be set by build.sh before sourcing:
#   - BUILD_EXT
#   - BUILD_DEBUG
#   - BUILD_MPI_REDUCTION

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
    CC=gcc CXX=g++

}
