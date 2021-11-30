#!/bin/sh
#$ -cwd
#$ -l {0[resource_group]}={0[node]}
#$ -l h_rt={0[time]}
#$ -N {0[JOB_NAME]}
#$ -e ../log/e
#$ -o ../log/o

# Load modules
. /etc/profile.d/modules.sh
module purge
module load intel
module load intel-mpi

echo start citylbm

export OMP_NUM_THREADS={0[NB_THREADS]}
export TMPDIR=/tmp

echo $PE_HOSTFILE
cp $PE_HOSTFILE hostfile.txt

ulimit -c unlimited

mpiexec.hydra -ppn {0[Processes_per_node]} \
       -n {0[mpi_procs]} \
       --bind-to board ../bin/citylbm.bdw_omp \
{0[CITYLBM_OPTIONS]}> ../log/logfile.txt 2>&1

echo end citylbm
