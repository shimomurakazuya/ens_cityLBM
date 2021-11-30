#!/bin/sh

#---------------------------------------------
# script used to submit citylbm job
# Usage ./job.sh
#---------------------------------------------
if ls ../bin/*.p100_cuda > /dev/null 2>&1; then
  qsub -g jh200053 batch_scripts/sub_p100_cuda.sh
elif ls ../bin/*.bdw_omp > /dev/null 2>&1; then
  qsub -g jh200053 batch_scripts/sub_bdw_omp.sh
elif ls ../bin/*.fugaku_a64fx_omp > /dev/null 2>&1; then
  pjsub batch_scripts/sub_fugaku_a64fx_omp.sh
elif ls ../bin/*.flow_a64fx_omp > /dev/null 2>&1; then
  pjsub batch_scripts/sub_flow_a64fx_omp.sh
else
  echo "No executable!"
fi
