#! /bin/bash

# module unload cray-hdf5-parallel
# module load cray-hdf5
module swap cpe cpe/23.03
module load cray-hdf5
module load cray-fftw
module load libfabric/1.20.1

if [ -d "${DEPSDIR}" ]; then
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DEPSDIR}/lib
fi
