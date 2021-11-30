#!/bin/sh
#$ -cwd
#$ -l f_node=9
#$ -l h_rt=16:00:00
#$ -N plant
#$ -v LD_LIBRARY_PATH=/apps/t3/sles12sp2/cuda/9.2.148/lib64:/apps/t3/sles12sp2/free/openmpi/2.1.2-opa10.9.0.1.2/gcc4.8.5/cuda9.2/lib/


rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe

. /etc/profile.d/modules.sh
module load cuda
module load openmpi

export OMP_NUM_THREADS=14
export TMPDIR=/tmp
export PSM2_MQ_RNDV_HFI_THRESH=128000

echo $PE_HOSTFILE
cp $PE_HOSTFILE hostfile.txt


mpirun -x PATH -x LD_LIBRARY_PATH -x PSM2_CUDA=1 -x PSM2_GPUDIRECT=1 -x CUDA_VISIBLE_DEVICES=0,1,2,3 -x HFI_UNIT=0,1,2,3 \
    -x PSM2_MULTIRAIL=2 \
    -mca orte_tmpdir_base /tmp \
    -mca coll_tuned_use_dynamic_rules 1 -mca coll_tuned_allgather_algorithm 2 \
    -mca btl_openib_cuda_rdma_limit 1000000 -mca mtl psm  --mca mpi_common_cuda_event_max 1000000 \
    -npernode 4 -n 36  --bind-to board  ./run.exe \
    -gpu_per_node              4 \
    -time_end                  24300.0 \
    -velocity_lbm              10.0  0.100 \
    -number_of_grid_point      80   2 \
    -domain_min                -2112.0  -2112.0    -8.0  \
    -domain_length              4224.0   4224.0  2560.0  \
    -cfr_steps                 750  750 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 99999 \
 > ../log/logfile.txt 2>&1
