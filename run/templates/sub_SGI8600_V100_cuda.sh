#!/bin/bash -ue
#PBS -q {0[resource_group]}
#PBS -l select={0[node]}:ncpus={0[nb_cores]}:mpiprocs={0[Processes_per_node]}:ompthreads={0[NB_THREADS]}:ngpus={0[Processes_per_node]}
#PBS -l walltime={0[time]}
#PBS -P CityLBM@{0[ACCOUNT]}
#PBS -N {0[JOB_NAME]}
#PBS -o ../log/
#PBS -j oe
echo start citylbm

if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi

export ofs_ensemble_idx={0[ofs_ensemble_idx]}

# check args
while getopts i: OPT
do
    case $OPT in 
        "i") export ofs_ensemble_idx=$OPTARG ;;
    esac
done

# env
. /etc/profile.d/modules.sh 
module purge 
module load gcc/7.4.0 mpt/2.23-ga cuda/11.0; export MPI_USE_CUDA=1

module list 2>&1 
export MPI_DSM_VERBOSE=1
export OMP_NUM_THREADS={0[NB_THREADS]}

env

#ulimit -c unlimited

echo start mpi

mpirun -np {0[mpi_procs]} omplace \
       ../bin/citylbm.v100_cuda \
{0[CITYLBM_OPTIONS]}|& tee ../log/logfile_$ofs_ensemble_idx.txt

echo finish mpi
