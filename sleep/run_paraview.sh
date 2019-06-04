#!/bin/sh
#$ -cwd
#$ -l f_node=1
#$ -l h_rt=0:05:00
#$ -N paraview_node

. /etc/profile.d/modules.sh
module load cuda openmpi paraview

cp $PE_HOSTFILE hostfile.txt

mpirun -np 4 --hostfile hostfile.txt pvserver --use-offscreen-rendering
