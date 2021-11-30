#!/bin/sh
#PJM -L "node={0[node]}"
#PJM -L "rscunit=fx"
#PJM -L "rscgrp={0[resource_group]}"
#PJM -L "elapse={0[time]}"
#PJM --mpi "proc={0[mpi_procs]}"
#PJM -s
#PJM -N {0[JOB_NAME]}

# Load modules
module purge
module load tcs/1.2.27

export OMP_NUM_THREADS={0[NB_THREADS]}
export OMPI_MCA_plm_ple_memory_allocation_policy=bind_local
export XOS_MMM_L_PAGING_POLICY=demand:demand:demand
export TMPDIR=/tmp

mkdir -p ../io/iofiles/restart ../io/output

mpiexec ../bin/citylbm.flow_a64fx_omp \
{0[CITYLBM_OPTIONS]}> ../log/logfile.txt 2>&1
