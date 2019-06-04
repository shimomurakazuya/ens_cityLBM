#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import math


def create_run(fname, arch, ncpus, gpu_per_node, number_of_grid_point):
    # parameters #
    num_nodes = int(math.ceil(ncpus/gpu_per_node)) # f_node
#    num_nodes = int(ncpus) # q_node
    
    ncpux = int(ncpus**0.5)
    ncpuy = int(ncpus/ncpux)

   

    # cal #
    velocity_lbm = [10.0, 0.05]
    if arch == 'gpu':
        cfr_steps = [ 500,  500 ]
    elif arch == 'cpu':
        cfr_steps = [ 25,  25 ]

    cfr_flags = [ 1,   1 ]
    
    restart_fs = [ 0, 99999 ] # [0] = 0 : false, 1 : true # [1] = restart step
    
    if arch == 'gpu':
        time_end0 = 20.1
        omp_num_threads = 14
        ncpus_cal = ncpus
        gpu_per_node_cal = gpu_per_node
    elif arch == 'cpu':
        time_end0 = 1.1
        omp_num_threads = 14
        ncpus_cal = ncpus
        gpu_per_node_cal = gpu_per_node
#        omp_num_threads = 28
#        ncpus_cal = int(ncpus/4)
#        gpu_per_node_cal = int(gpu_per_node/4)

    domain_min0    = [  0.0,  0.0,  -4.0 ]
    domain_length0 = [  288.0,   288.0,  896.0 ]
    
    
    time_end      = time_end0
    domain_min    = [ domain_min0   [0]*ncpux, domain_min0   [1]*ncpuy, domain_min0   [2] ]
    domain_length = [ domain_length0[0]*ncpux, domain_length0[1]*ncpuy, domain_length0[2] ]


    # strings
    strings  = '#!/bin/sh\n'
    strings += '\n'

    strings += '#$ -cwd\n'
    strings += '#$ -l f_node=' + str(num_nodes) + '\n'
#    strings += '#$ -l q_node=' + str(num_nodes) + '\n'
    strings += '#$ -l h_rt=0:30:00\n'
    strings += '#$ -N ws_' + arch + str(ncpus) + '\n'
    strings += '#$ -p -4\n'
    strings += '#$ -v LD_LIBRARY_PATH=/apps/t3/sles12sp2/free/openmpi/1.10.7/gnu-cuda8.0/lib:/apps/t3/sles12sp2/cuda/8.0.61/lib.375.66:/apps/t3/sles12sp2/cuda/8.0.61/lib64:/apps/t3/sles12sp2/free/hdf5/1.10.1/gnu/lib\n'
    strings += '\n'


    strings += 'rm run.exe\n'
    strings += 'ln -s ../bin/run.exe\n'
    strings += 'chmod 777 run.exe\n'
    strings += '\n'


    strings += '. /etc/profile.d/modules.sh\n'
    strings += 'module load cuda/8.0.61\n'
    strings += 'module load /apps/t3/sles12sp2/modules/modulefiles/gpu-direct/openmpi/1.10.7/gnu-cuda8.0\n'
    strings += '\n'


    strings += 'export OMP_NUM_THREADS=' + str(omp_num_threads) +'\n'
    strings += '\n'

    strings += 'cp $PE_HOSTFILE hostfile.txt\n'
    strings += '\n'

    # mpirun
    strings += 'mpirun -x PATH -x LD_LIBRARY_PATH -x PSM2_CUDA=1 -x PSM2_GPUDIRECT=1 -x CUDA_VISIBLE_DEVICES=0,1,2,3 -x HFI_UNIT=1'  + ' \\\n'
    strings += '    -mca orte_tmpdir_base /tmp'  + ' \\\n'
    strings += '    -mca coll_tuned_use_dynamic_rules 1 -mca coll_tuned_allgather_algorithm 2'  + ' \\\n'
    strings += '    -mca btl_openib_cuda_rdma_limit 1000000 -mca mtl psm  --mca mpi_common_cuda_event_max 1000000'  + ' \\\n'
    strings += '    -npernode ' + str(gpu_per_node) + ' -n ' + str(ncpus_cal) + '  --bind-to board   ./run.exe' +   ' \\\n'
    strings += '    -gpu_per_node              ' + str(gpu_per_node_cal) + ' \\\n'
    strings += '    -time_end                  ' + str(time_end) + ' \\\n'
    strings += '    -velocity_lbm              ' + str(velocity_lbm[0]) + ' ' + str(velocity_lbm[1]) + ' \\\n'
    strings += '    -number_of_grid_point      ' + str(number_of_grid_point[0]) + ' ' + str(number_of_grid_point[1]) + ' \\\n'
    strings += '    -domain_min                ' + str(domain_min[0]) + ' ' + str(domain_min[1]) + ' ' + str(domain_min[2]) + ' \\\n'
    strings += '    -domain_length             ' + str(domain_length[0]) + ' ' + str(domain_length[1]) + ' ' + str(domain_length[2]) + ' \\\n'
    strings += '    -cfr_steps                 ' + str(cfr_steps[0]) + ' ' + str(cfr_steps[1]) + ' \\\n'
    strings += '    -cfr_flags                 ' + str(cfr_flags[0]) + ' ' + str(cfr_flags[1]) + ' \\\n'
    strings += '    -restart_flags_and_step    ' + str(restart_fs[0]) + ' ' + str(restart_fs[1]) + ' \\\n'
    strings += ' > ' + '../log/logfile.txt 2>&1'


    # fout #
    f = open(fname, 'w')
    f.write(strings)
    f.close()



if __name__ == '__main__':
    rundir = os.getcwd()
    print(rundir)
    
    fname_run = 'run_amr_performance.sh'


    # argument
    parser = argparse.ArgumentParser()
    parser.add_argument("ncpus", help="number of processes")
    parser.add_argument("arch",  help="cpu or gpu")
    args = parser.parse_args()


    # parameters
    ncpus        = int(args.ncpus)
    arch         = str(args.arch)

    gpu_per_node = 4

    number_of_grid_point = [ 56, 2 ] # leaf resolution  #

    # python function
    create_run(fname_run, arch, ncpus, gpu_per_node, number_of_grid_point)


    # shell
    os.chmod(fname_run, 0o777)

