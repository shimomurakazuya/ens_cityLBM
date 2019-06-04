#!/bin/sh
#$ -cwd
#$ -l q_node=4
#$ -l h_rt=0:05:00
#$ -N test_lbm
#$ -p -3
#$ -v LD_LIBRARY_PATH=/apps/t3/sles12sp2/free/openmpi/2.1.2/gnu/lib:/apps/t3/sles12sp2/cuda/8.0.61/lib.375.66:/apps/t3/sles12sp2/cuda/8.0.61/lib64:/apps/t3/sles12sp2/free/hdf5/1.10.1/gnu/lib


rm run.exe
ln -s ../bin/run.exe
chmod 777 run.exe

. /etc/profile.d/modules.sh
module load cuda/8.0.61
module load openmpi/2.1.2

export OMP_NUM_THREADS=7

echo $PE_HOSTFILE
cp $PE_HOSTFILE hostfile.txt



#mpirun  -npernode 1 -n 4 nvprof  --annotate-mpi openmpi  --profile-all-processes   --events all --metrics all  -o simpleMPI.%q{OMPI_COMM_WORLD_RANK}.nvprof  ./run.exe \
#mpirun  -npernode 1 -n 4 nvprof  --events all --metrics all  -o simpleMPI2.%q{OMPI_COMM_WORLD_RANK}.nvprof  ./run.exe \
#mpirun  -npernode 1 -n 4 nvprof  --events all --metrics all  --unified-memory-profiling per-process-device  -o simpleMPI.%q{OMPI_COMM_WORLD_RANK}.nvprof  ./run.exe \
#mpirun  -npernode 4 -n 4 nvprof  --events all --metrics all  --unified-memory-profiling per-process-device  -o simpleMPI.%q{OMPI_COMM_WORLD_RANK}.nvprof  ./run.exe \

#mpirun  -npernode 4 -n 4 nvprof   --unified-memory-profiling per-process-device --metrics all  --events all  ./run.exe \
#mpirun  -npernode 4 -n 4 nvprof  --cpu-profiling on  -f -o simpleMPI.%q{OMPI_COMM_WORLD_RANK}.nvprof  ./run.exe \
#mpirun  -npernode 4 -n 4 nvprof --cpu-profiling on --cpu-thread-tracing on  --unified-memory-profiling per-process-device --metrics all   --events all -f -o simpleMPI.%q{OMPI_COMM_WORLD_RANK}.nvprof  ./run.exe \
#mpirun  -npernode 4 -n 4 nvprof --metrics all   --events all -f -o simpleMPI.%q{OMPI_COMM_WORLD_RANK}.nvprof  ./run.exe \

# simple timeline #
#mpirun  -npernode 1 -n 4  nvprof -f -o simpleMPI.%q{OMPI_COMM_WORLD_RANK}.nvprof ./run.exe \

# cal : ok #
#mpirun  -npernode 1 -n 4  ./run.exe \

# flops : ok #
#mpirun  -npernode 1 -n 4 nvprof   --unified-memory-profiling per-process-device --metrics all  --events all  ./run.exe \


# detail of kernel function #
#mpirun  -npernode 1 -n 4  -x LD_LIBRARY_PATH -x PATH   nvprof  --events all --metrics all -f -o simpleMPI.%q{OMPI_COMM_WORLD_RANK}.nvprof ./run.exe \

mpirun  -npernode 1 -n 4  -x LD_LIBRARY_PATH -x PATH   nvprof  -f -o simpleMPI.%q{OMPI_COMM_WORLD_RANK}.nvprof ./run.exe \
    -gpu_per_node              4 \
    -time_end                  0.004 \
    -velocity_lbm              2.0  0.05 \
    -number_of_grid_point      12 2 \
    -domain_min                -4.8  -1.2 -0.2  \
    -domain_length              9.6   2.4  2.4  \
    -cfr_steps                 50  50 \
    -cfr_flags                 1  1 \
    -restart_flags_and_step    0 99999 \
 > ../log/logfile.txt 2>&1

