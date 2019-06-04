#!/bin/sh

#$ -cwd
#$ -l q_node=1
#$ -l h_rt=0:15:00
#$ -N test_lbm
#$ -p -4
#$ -v LD_LIBRARY_PATH=/apps/t3/sles12sp2/free/openmpi/2.1.2-cuda-hfi/gnu/lib:/apps/t3/sles12sp2/cuda/8.0.61/lib.375.66:/apps/t3/sles12sp2/cuda/8.0.61/lib64:/apps/t3/sles12sp2/free/hdf5/1.10.1/gnu/lib

rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe

. /etc/profile.d/modules.sh
module load cuda/8.0.61
module load openmpi/2.1.2

export OMP_NUM_THREADS=14
export PSM2_MQ_RNDV_HFI_THRESH=128000

cp $PE_HOSTFILE hostfile.txt

mpirun -x PATH -x LD_LIBRARY_PATH -x PSM2_CUDA=1 -x PSM2_GPUDIRECT=1 -x CUDA_VISIBLE_DEVICES=0,1,2,3 -x HFI_UNIT=0,1,2,3 \
    -x PSM2_MULTIRAIL=2 \
    -mca orte_tmpdir_base /tmp \
    -mca coll_tuned_use_dynamic_rules 1 -mca coll_tuned_allgather_algorithm 2 \
    -mca btl_openib_cuda_rdma_limit 1000000 -mca mtl psm  --mca mpi_common_cuda_event_max 1000000 \
    -npernode 1 -n 1   --bind-to board  ./run.exe \
    -gpu_per_node              4 \
    -time_end                  600.0 \
    -velocity_lbm              1.5  0.1 \
    -number_of_grid_point      36  0 \
    -domain_min                -0.51  -0.51  -0.00001 \
    -domain_length              1.02   1.02   0.00002 \
    -cfr_steps                 5000  25000 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 0 \
 > ../log/logfile.txt 2>&1
