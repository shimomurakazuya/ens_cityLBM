#!/bin/bash
#SBATCH -J {0[JOB_NAME]}
#SBATCH -p {0[resource_group]}
#SBATCH -w {0[nodelist]}
#SBATCH -N {0[node]}
#SBATCH -n {0[mpi_procs]}
#SBATCH --cpus-per-task={0[nb_cores]}
#SBATCH --cpus-per-gpu={0[nb_cores]}
#SBATCH --gpus-per-node={0[Processes_per_node]}
#SBATCH --time {0[time]}
#SBATCH -o ./stdout.%J
#SBATCH -e ./stderr.%J

# Load modules
# default: module load gcc/8.3.1 cuda/11.2 ompi-cuda/4.1.1-11.2
module purge
{0[MODULE_SETTINGS]}

# Python environment (optional)
{0[PYTHON_SETTINGS]}

# Environmental variables
# default: export OMP_NUM_THREADS=32
{0[ENV_SETTINGS]}

# Execute app
{0[EXEC_COMMANDS]}
