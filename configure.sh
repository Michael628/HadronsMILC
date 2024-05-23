#! /bin/bash

# Build the hadrons app

SRCDIR="."
HADRONSINSTALLDIR="../install-hadrons-avx512-knl"
GRIDINSTALLDIR="../install-grid-avx512-knl"
PK_CXX="CC"
PK_CC="cc"

./bootstrap.sh

${SRCDIR}/configure \
  --with-hadrons=${HADRONSINSTALLDIR}      \
  --with-grid=${GRIDINSTALLDIR}  \
  --host=x86_64-unknown-linux-gnu \
  CXX=${PK_CXX}              \
  CXXFLAGS="-O3"



