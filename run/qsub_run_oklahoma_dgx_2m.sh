#!/bin/sh
#PBS -q full
#PBS -l select=1:ncpus=96:mpiprocs=16:ompthreads=6
#PBS -P test
#PBS -l walltime=16:30:00
#PBS -N plant

cd $PBS_O_WORKDIR

. /etc/profile.d/modules.sh
module load cuda/10.0
module load compiler/gcc-7.3.0
module load mvapich2/2.3/gcc-7.3.0.lp


rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe


export MV2_SMP_USE_CMA=0
#export MV2_GPUDIRECT_GDRCOPY_LIB=/home/app/mvapich2/2.2-gdr-cuda8.0/gnu/gdrcopy-master/libgdrapi.so
export MV2_GPUDIRECT_GDRCOPY_LIB=/opt/mvapich2/2.3/gdrcopy-1.3/libgdrapi.so
export MV2_USE_CUDA=1
export MV2_USE_GPUDIRECT_RDMA=1
export MV2_GPUDIRECT_GDRCOPY=1
export MV2_USE_GPUDIRECT_GDRCOPY=1
export MV2_CUDA_IPC=1
export MV2_CUDA_ENABLE_MANAGED=1   # this flag must be set of managed memory is used otherwise it fails
export MV2_CUDA_MANAGED_IPC=1

export OMP_NUM_THREADS=3


env > ../log/logfile.txt 2>&1

mpiexec  \
    -n 16  --bind-to board   ./run.exe \
    -gpu_per_node              16 \
    -time_end                  24300.0 \
    -velocity_lbm              10.0  0.100 \
    -number_of_grid_point      80   2 \
    -domain_min                -2048.0  -2048.0    -8.0  \
    -domain_length              4096.0   4096.0  2560.0  \
    -cfr_steps                 750  750 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 99999 \
 >> ../log/logfile.txt 2>&1
