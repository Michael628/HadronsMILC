module purge

module load gompi ucc_cuda fftw hdf5

if [ -d "${DEPSDIR}" ]; then
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DEPSDIR}/lib
fi
