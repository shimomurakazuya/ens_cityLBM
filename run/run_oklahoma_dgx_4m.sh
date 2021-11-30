#!/bin/sh

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


mpiexec  \
    -n 16  --bind-to board   ./run.exe \
    -gpu_per_node              16 \
    -time_end                  23400.0 \
    -velocity_lbm              10.0  0.075 \
    -number_of_grid_point      40   2 \
    -domain_min                -3840.0  -3840.0   -16.0  \
    -domain_length              7680.0   7680.0  2560.0  \
    -cfr_steps                 500  500 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 99999 \
 > ../log/logfile.txt 2>&1


