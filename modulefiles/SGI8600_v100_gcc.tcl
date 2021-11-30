#%Module

## setup intel compiler
module load gcc/7.4.0 
 
## setup cuda and mpi
module load cuda/11.0
module load mpt/2.23-ga
setenv MPI_USE_CUDA 1

## setup CityLBM::DEVICE
setenv DEVICE "v100"
