#!/bin/sh
#$ -cwd
#$ -l f_node=16
#$ -l h_rt=00:02:00
#$ -N test_lbm
#$ -p -5
#$ -v LD_LIBRARY_PATH=/apps/t3/sles12sp2/free/openmpi/1.10.7/gnu-cuda8.0/lib:/apps/t3/sles12sp2/cuda/8.0.61/lib.375.66:/apps/t3/sles12sp2/cuda/8.0.61/lib64:/apps/t3/sles12sp2/free/hdf5/1.10.1/gnu/lib


rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe

. /etc/profile.d/modules.sh
module load cuda/8.0.61
module load /apps/t3/sles12sp2/modules/modulefiles/gpu-direct/openmpi/1.10.7/gnu-cuda8.0

export OMP_NUM_THREADS=14
export TMPDIR=/tmp
export PSM2_MQ_RNDV_HFI_THRESH=128000

echo $PE_HOSTFILE
cp $PE_HOSTFILE hostfile.txt


mpirun -x PATH -x LD_LIBRARY_PATH -x PSM2_CUDA=1 -x PSM2_GPUDIRECT=1 -x CUDA_VISIBLE_DEVICES=0,1,2,3 -x HFI_UNIT=0 \
    -mca orte_tmpdir_base /tmp \
    -mca coll_tuned_use_dynamic_rules 1 -mca coll_tuned_allgather_algorithm 2 \
    -mca btl_openib_cuda_rdma_limit 1000000 -mca mtl psm  --mca mpi_common_cuda_event_max 1000000 \
    -npernode 4 -n 64  --bind-to board  ./run.exe \
    -gpu_per_node              4 \
    -time_end                  20.1 \
    -velocity_lbm              10.0  0.05 \
    -number_of_grid_point      56   2 \
    -domain_min                   0     0   -4  \
    -domain_length             576  576  896  \
    -cfr_steps                 500  5000 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 99999 \
 > ../log/logfile.txt 2>&1


