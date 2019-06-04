#!/bin/sh

rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe

qsub \
  -q gp28 \
  -P lbm \
  -l select=1:ncpus=28:mpiprocs=4:ompthreads=7 \
  -l walltime=00:02:00 \
  -N test \
run_amr.sh
