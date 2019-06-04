#!/bin/sh
#$ -cwd
#$ -l q_node=4
#$ -l h_rt=00:05:00
#$ -N test
#$ -p -5
#$ -v LD_LIBRARY_PATH=/apps/t3/sles12sp2/free/openmpi/1.10.7/gnu-cuda8.0/lib:/apps/t3/sles12sp2/cuda/8.0.61/lib.375.66:/apps/t3/sles12sp2/cuda/8.0.61/lib64:/apps/t3/sles12sp2/free/hdf5/1.10.1/gnu/lib


rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe

. /etc/profile.d/modules.sh
module load cuda/8.0.61
module load /apps/t3/sles12sp2/modules/modulefiles/gpu-direct/openmpi/1.10.7/gnu-cuda8.0

export OMP_NUM_THREADS=7

echo $PE_HOSTFILE
cp $PE_HOSTFILE hostfile.txt


mpirun -x PATH -x LD_LIBRARY_PATH -x PSM2_CUDA=1 -x PSM2_GPUDIRECT=1 -x CUDA_VISIBLE_DEVICES=0 -x HFI_UNIT=1 \
    --mca btl_openib_cuda_rdma_limit 1000000 -mca mtl psm  --mca mpi_common_cuda_event_max 1000000 \
    -npernode 1 -n 4   --bind-to board  ./run.exe \
    -gpu_per_node              4 \
    -time_end                  660.0 \
    -velocity_lbm              15.0  0.10 \
    -number_of_grid_point      10   2 \
    -domain_min                -6.6  -3.3  -1.1  \
    -domain_length             13.2   6.6   2.2    \
    -cfr_steps                 1000  100000 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 0 \
 > ../log/logfile.txt 2>&1



