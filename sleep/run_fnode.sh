#!/bin/sh
#$ -cwd
#$ -l f_node=1
#$ -l h_rt=24:00:00
#$ -N ond_node

. /etc/profile.d/modules.sh
module load intel
module load cuda
module load openmpi

cp $PE_HOSTFILE hostfile.txt

sleep 864000
