#!/bin/bash
#PJM -L "node={0[node]}"
#PJM -L "rscgrp={0[resource_group]}"
#PJM -L "elapse={0[time]}"
#PJM -s
#PJM -g {0[ACCOUNT]}
#PJM --mpi proc={0[mpi_procs]}

# Load modules
# default: module load gcc/8.3.1 cuda/11.2 ompi-cuda/4.1.1-11.2
module purge
{0[LOAD_MODULES]}

# environmental variables
# default: export UCX_MEMTYPE_CACHE=n; export UCX_IB_GPU_DIRECT_RDMA=no
{0[ENV_VARIABLES]}

echo start citylbm
mpiexec -machinefile $PJM_O_NODEINF -np $PJM_MPI_PROC -npernode {0[Processes_per_node]} ../bin/citylbm.A100_cuda \
{0[CITYLBM_OPTIONS]}> ../log/logfile.txt 2>&1

echo end citylbm
