#!/bin/sh
#$ -cwd
#$ -l rt_F=9
#$ -l h_rt=00:15:00
#$ -N oklahoma4m


rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe

source /etc/profile.d/modules.sh
module load cuda/9.2/9.2.88.1
module load openmpi

export OMP_NUM_THREADS=14
export TMPDIR=/tmp
export PSM2_MQ_RNDV_HFI_THRESH=128000

echo $PE_HOSTFILE
cp $PE_HOSTFILE hostfile.txt


mpirun  \
    -x PSM2_CUDA=1 -x PSM2_GPUDIRECT=1 -x CUDA_VISIBLE_DEVICES=0,1,2,3\
    -mca orte_tmpdir_base /tmp \
    -mca coll_tuned_use_dynamic_rules 1 -mca coll_tuned_allgather_algorithm 2 \
    -mca btl_openib_cuda_rdma_limit 1000000 -mca mtl psm  --mca mpi_common_cuda_event_max 1000000 \
    -npernode 4 -n 36  --bind-to board  ./run.exe \
    -gpu_per_node              4 \
    -time_end                  36000.0 \
    -velocity_lbm              10.0  0.075 \
    -number_of_grid_point      40   2 \
    -domain_min                -3840.0  -3840.0   -16.0  \
    -domain_length              7680.0   7680.0  2560.0  \
    -cfr_steps                 2000  2000 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 99999 \
 > ../log/logfile.txt 2>&1

