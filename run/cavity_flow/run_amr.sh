#!/bin/sh

cd $PBS_O_WORKDIR

source /home/usr6/12IH0099/env/env.sh

mpirun  --bind-to board  -np 4 -hostfile $PBS_NODEFILE  -x OMP_NUM_THREADS=$OMP_NUM_THREADS  ./run.exe \
    -time_end                  50.0 \
    -velocity_lbm              1.0 0.15 \
    -number_of_grid_point      8 0 \
    -domain_min                -0.52 -0.52 -0.052 \
    -domain_length             1.04 1.04 0.12 \
    -cfr_steps                 100 100 10000 \
    -cfr_flags                 1 1 0 \
    -restart_flags_and_step    0 99999 \
 > ../log/logfile.txt 2>&1
