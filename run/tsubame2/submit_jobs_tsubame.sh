#!/bin/sh

rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe

t2sub \
  -N test_lbm \
  -q S \
  -r n \
  -p 2 \
  -et 0 \
  -W group_list=t2g-jh170031 \
  -l walltime=00:01:00 \
  -l select=4:ncpus=24:mpiprocs=1:ompthreads=24:mem=18gb:gpus=3 \
run_amr.sh
