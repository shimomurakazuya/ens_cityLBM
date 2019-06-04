#!/bin/sh
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=0:07:00
#$ -N test_lbm
#$ -p -3
#$ -v LD_LIBRARY_PATH=/apps/t3/sles12sp2/free/openmpi/1.10.7/gnu-cuda8.0/lib:/apps/t3/sles12sp2/cuda/8.0.61/lib.375.66:/apps/t3/sles12sp2/cuda/8.0.61/lib64:/apps/t3/sles12sp2/free/hdf5/1.10.1/gnu/lib


rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe

. /etc/profile.d/modules.sh
module load cuda/8.0.61
module load /apps/t3/sles12sp2/modules/modulefiles/gpu-direct/openmpi/1.10.7/gnu-cuda8.0

export OMP_NUM_THREADS=1

echo $PE_HOSTFILE
cp $PE_HOSTFILE hostfile.txt

# ok : flops #
mpirun  -npernode 1 -n 4 nvprof   --unified-memory-profiling per-process-device --metrics all  --events all  ./run.exe \
    -gpu_per_node              4 \
    -time_end                  0.004 \
    -velocity_lbm              2.0  0.05 \
    -number_of_grid_point      6 2 \
    -domain_min                -4.8  -1.2 -0.2  \
    -domain_length              9.6   2.4  2.4  \
    -cfr_steps                 50  50 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 99999 \
 > ../log/logfile.txt 2>&1


