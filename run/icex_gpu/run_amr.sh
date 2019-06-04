#!/bin/sh

cd $PBS_O_WORKDIR

. /etc/profile.d/modules.sh
module load cuda/8.0 intel/17.0.2.174

source /home/g5/a170085/env/env.sh

export MV2_USE_CUDA=1


mpiexec ./run.exe \
    -gpu_per_node              4 \
    -time_end                  50.0 \
    -velocity_lbm              1.0 0.15 \
    -number_of_grid_point      8 0 \
    -domain_min                -0.52 -0.52 -0.052 \
    -domain_length             1.04 1.04 0.12 \
    -cfr_steps                 100 200 \
    -cfr_flags                 1 1 \
    -restart_flags_and_step    0 99999 \
 > ../log/logfile.txt 2>&1
