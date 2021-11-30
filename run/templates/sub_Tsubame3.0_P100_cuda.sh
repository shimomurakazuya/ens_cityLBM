#!/bin/sh
#$ -cwd
#$ -l {0[resource_group]}={0[node]}
#$ -l h_rt={0[time]}
#$ -N {0[JOB_NAME]}
#$ -e ../log/e
#$ -o ../log/o

# Load modules
# default: module load cuda/11.2.146 gcc/10.2.0-cuda openmpi/3.1.4-opa10.10-t3
. /etc/profile.d/modules.sh
module purge
{0[LOAD_MODULES]}

echo start citylbm

export OMP_NUM_THREADS={0[NB_THREADS]}

# environmental variables
# default export PSM2_MQ_RNDV_HFI_THRESH=128000
export TMPDIR=/tmp
{0[ENV_VARIABLES]}

echo $PE_HOSTFILE
cp $PE_HOSTFILE hostfile.txt

ulimit -c unlimited

mpirun -x PATH -x LD_LIBRARY_PATH -x PSM2_CUDA=1 -x PSM2_GPUDIRECT=1 -x CUDA_VISIBLE_DEVICES=0,1,2,3 -x HFI_UNIT=0,1,2,3 \
       -x PSM2_MULTIRAIL=2 \
       -mca orte_tmpdir_base /tmp \
       -mca coll_tuned_use_dynamic_rules 1 -mca coll_tuned_allgather_algorithm 2 \
       -mca btl_openib_cuda_rdma_limit 1000000 -mca mtl psm  --mca mpi_common_cuda_event_max 1000000 \
       -npernode {0[Processes_per_node]} \
       -n {0[mpi_procs]} \
       --bind-to board ../bin/citylbm.p100_cuda \
{0[CITYLBM_OPTIONS]}> ../log/logfile.txt 2>&1

echo end citylbm
