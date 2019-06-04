#!/bin/sh
#PBS -q l-small
#PBS -l select=1:mpiprocs=4:ompthreads=9
#PBS -W group_list=gi37
#PBS -l walltime=00:04:00
cd $PBS_O_WORKDIR

. /etc/profile.d/modules.sh

module load cuda9/9.0.176
module load mvapich2/gdr/2.3a/gnu 

export MV2_USE_CUDA=1 # device
export MV2_USE_GPUDIRECT=1 # device
export MV2_GPUDIRECT_GDRCOPY_LIB=/lustre/app/mvapich2-gdr/ofed4.2/cuda9.0/gdrcopy/gnu/lib64/libgdrapi.so
export MV2_CUDA_IPC=1 # if this flag is 0 NVLINK is not used
#export MV2_CUDA_SMP_IPC=1
#export MV2_CUDA_BLOCK_SIZE=131072
#export MV2_CUDA_BLOCK_SIZE=262144 # gives best results so far
#export MV2_CUDA_BLOCK_SIZE=393216
#export MV2_CUDA_BLOCK_SIZE=524288
#export MV2_GPUDIRECT_RECEIVE_LIMIT=16384
#export MV2_GPUDIRECT_LIMIT=262144
#export MV2_USE_SMP_GDR=1
export MV2_ENABLE_AFFINITY=0
#export MV2_CUDA_IPC_THRESHOLD=262144
#export MV2_GPUDIRECT_LIMIT=4194304

export MV2_CUDA_ENABLE_MANAGED=1 # managed memory
export MV2_CUDA_MANAGED_IPC=1 # managed memory

#export MV2_USE_MCAST=1  # not supported
#export MV2_MCAST_NUM_NODES_THRESHOLD=16 # not supported

#export KMP_AFFINITY=verbose,granularity=fine,scatter,0,0
#export I_MPI_PIN_DOMAIN=omp:compact
#export I_MPI_PIN_CELL=core

OUTFILE="output.log"
#NPROF+="--aggregate-mode off --event-collection-mode continuous"
#NPROF+=" --metrics nvlink_total_data_transmitted"
#NPROF+=" --metrics nvlink_total_data_received"
#NPROF+=" --metrics nvlink_transmit_throughput"
#NPROF+=" --metrics nvlink_receive_throughput"

NPROF+="--metrics flop_count_dp "
NPROF+="--metrics dram_read_transactions --metrics dram_write_transactions "
#NPROF+="--metrics gld_transactions       --metrics gst_transactions "

#mpiexec ./a.out 12 32 > $OUTFILE 2>&1


mpiexec -n 4 \
    --bind-to board  ../bin/run.exe \
    -gpu_per_node              4 \
    -time_end                  30.1 \
    -velocity_lbm              10.0  0.05 \
    -number_of_grid_point      56   2 \
    -domain_min                   0     0   -4  \
    -domain_length             640  640  896  \
    -cfr_steps                 500  5000 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 99999 \
 > ../log/logfile.txt 2>&1


