#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import math


base_ncpux = 4
    

def create_run(fname, arch, ncpus, gpu_per_node, number_of_grid_point):
    # parameters #
    num_nodes = int(math.ceil(ncpus/gpu_per_node)) # f_node
#    num_nodes = int(ncpus) # q_node

    queue_type = 'l-small'
    if num_nodes >= 9:
        queue_type = 'l-large'
    elif num_nodes >= 5:
        queue_type = 'l-medium'
    else:
        queue_type = 'l-small'
    
    ncpux = int(ncpus**0.5)
    ncpuy = int(ncpus/ncpux)

    MV2_USE_GPUDIRECT01 = 0
    if num_nodes == 1:
        MV2_USE_GPUDIRECT01 = 1
   

    # cal #
    velocity_lbm = [10.0, 0.05]
    if arch == 'gpu':
        cfr_steps = [ 500,  500 ]
    elif arch == 'cpu':
        cfr_steps = [ 25,  25 ]

    cfr_flags = [ 1,   1 ]
    
    restart_fs = [ 0, 99999 ] # [0] = 0 : false, 1 : true # [1] = restart step
    
    if arch == 'gpu':
        time_end0 = 30.1
        ncpus_cal = ncpus
        gpu_per_node_cal = gpu_per_node
    elif arch == 'cpu':
        time_end0 = 1.1
        ncpus_cal = ncpus
        gpu_per_node_cal = gpu_per_node
#        ncpus_cal = int(ncpus/4)
#        gpu_per_node_cal = int(gpu_per_node/4)

    domain_min0    = [  0.0,  0.0,  -4.0 ]
#    domain_length0 = [  288.0,   288.0,  896.0 ]
    domain_length0 = [  320.0,   320.0,  896.0 ]
#    domain_length0 = [  336.0,   336.0,  896.0 ]
    
    time_end      = time_end0
    # strong 
    domain_min    = [ domain_min0   [0]*base_ncpux, domain_min0   [1]*base_ncpux, domain_min0   [2] ]
    domain_length = [ domain_length0[0]*base_ncpux, domain_length0[1]*base_ncpux, domain_length0[2] ]


    # strings
    strings  = '#!/bin/sh\n'
    strings += '#PBS -q ' + queue_type + '\n'
    strings += '#PBS -l select=' + str(num_nodes) + ':mpiprocs=4:ompthreads=9\n'
    strings += '#PBS -W group_list=gi37\n'
    strings += '#PBS -l walltime=00:15:00\n'
    strings += 'cd $PBS_O_WORKDIR\n'
    strings += '\n'

    strings += '. /etc/profile.d/modules.sh\n'
    strings += 'module load cuda9/9.0.176\n'
    strings += 'module load mvapich2/gdr/2.3a/gnu\n'
    strings += '\n'

    strings += 'export MV2_USE_CUDA=1\n'
    strings += 'export MV2_USE_GPUDIRECT=' + str(MV2_USE_GPUDIRECT01) + '\n'
    strings += 'export MV2_GPUDIRECT_GDRCOPY_LIB=/lustre/app/mvapich2-gdr/ofed4.2/cuda9.0/gdrcopy/gnu/lib64/libgdrapi.so\n'
    strings += 'export MV2_CUDA_IPC=1\n'
    strings += 'export MV2_ENABLE_AFFINITY=0\n'
    strings += 'export MV2_CUDA_ENABLE_MANAGED=1\n'
    strings += 'export export MV2_CUDA_MANAGED_IPC=1\n'
    strings += '\n'


    # mpirun
    strings += 'mpiexec -n ' + str(ncpus) + ' --bind-to board  ../bin/run.exe \\\n'
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



