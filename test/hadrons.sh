#! /bin/bash

# Set run environment
eval $(pyfm workspace env --system scalar)
if [ $? -ne 0 ]; then
  echo "Error: failed to source system env for 'PYFM_SYSTEM_NAME'" >&2
  exit 1
fi

echo $PATH
# PYFM_EXECUTABLE selects the application; one of:
#   grid_lma | HadronsMILC | make_links_hisq
# Defaults to HadronsMILC if not set in the environment.
executable=${PYFM_EXECUTABLE:-HadronsMILC}

export OPT="${OPT} --comms-overlap"
runargs=" --grid $LATTICE --mpi ${LAYOUT} $OPT"

for input in ${INPUTLIST}
do
  output=out/$(basename "${input%.xml}").${SLURM_JOB_ID}
  echo "Input file ${input}"
  echo "Output file ${output}"
  echo "START_RUN `date`" >> ${output}
  echo "RUNARGS = ${runargs}" >> ${output}

  date >> ${output}
  echo "Executable: $(which ${executable})" >> ${output}
  argstr="${input} ${runargs}"
  export APP="${executable} ${argstr}"
  echo ${APP} >> ${output}
  cmd="${APP}"

  echo ${cmd} >> ${output}
  ${cmd} >> ${output} &
  # ${cmd}
  echo ${cmd}
done

wait
exit 0
