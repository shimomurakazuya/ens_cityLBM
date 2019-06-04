#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

#home = '/home/g5/a170085'
home = '/home/usr6/12IH0099' # tsubame

rundir = os.getcwd()
#bindir = rundir + '/../bin'
#logdir = rundir + '/../log'

print(rundir)

fname_pbs = 'submit_jobs_tsubame.sh'
fname_run = 'run_amr.sh'

# parameters #
# pbs #
ncpus = 24
num_nodes = 9
mpiprocs_node = 1


# cal #
time_end = 50.0
velocity_lbm = [1.0, 0.15]
number_of_grid_point = [ 9, 0 ] # leaf resolution  #
cfr_steps = [ 100,   200,   10000 ]
cfr_flags = [ 1,   1,   0]

restart_fs = [ 0, 99999 ] # [0] = 0 : false, 1 : true # [1] = restart step

domain_min    = [ -0.52,  -0.52,  -0.01375 ]
domain_length = [  1.04,   1.04,   0.02750 ]

dx_leaf = domain_length[number_of_grid_point[1]] / number_of_grid_point[0]

domain_min[2]    = -dx_leaf*0.5
domain_length[2] =  dx_leaf


def create_pbs_script(fname):
    strings  = '#!/bin/sh\n'
    strings += '\n'

    strings += 'rm run.exe\n'
    strings += 'ln -s ../bin/run.exe\n'
    strings += 'chmod 777 run.exe\n'
    strings += '\n'

    strings += 't2sub \\\n'
    strings += '  -N test_lbm \\\n'
    strings += '  -q S \\\n'
    strings += '  -r n \\\n'
    strings += '  -p 0 \\\n'
    strings += '  -et 1 \\\n'
    strings += '  -W group_list=t2g-jh170031 \\\n'
    strings += '  -l walltime=00:02:00 \\\n'
#    strings += '  -l select=' + str(num_nodes) + ':ncpus=' + str(ncpus) + ':mpiprocs=' + str(mpiprocs_node) + ':ompthreads=' + str(int(ncpus/mpiprocs_node)) + ' \\\n'
    strings += '  -l select=' + str(num_nodes) + ':ncpus=' + str(ncpus) + ':mpiprocs=' + str(mpiprocs_node) + ':ompthreads=' + str(int(ncpus/mpiprocs_node)) + ':mem=18gb:gpus=3 \\\n'
    strings + '\n'

    strings += fname_run

    # fout #
    f = open(fname, 'w')
    f.write(strings)
    f.close()


def create_run(fname):
    strings  = '#!/bin/sh\n'
    strings += '\n'

    strings += 'cd $PBS_O_WORKDIR\n'
    strings += '\n'
    strings += 'source ' + home + '/env/env.sh\n'
    strings += '\n'


    strings += 'mpirun  --bind-to board  -np ' + str(int(num_nodes*mpiprocs_node)) + ' -hostfile $PBS_NODEFILE  -x OMP_NUM_THREADS=$OMP_NUM_THREADS  ./run.exe' +   ' \\\n'
    strings += '    -time_end                  ' + str(time_end) + ' \\\n'
    strings += '    -velocity_lbm              ' + str(velocity_lbm[0]) + ' ' + str(velocity_lbm[1]) + ' \\\n'
    strings += '    -number_of_grid_point      ' + str(number_of_grid_point[0]) + ' ' + str(number_of_grid_point[1]) + ' \\\n'
    strings += '    -domain_min                ' + str(domain_min[0]) + ' ' + str(domain_min[1]) + ' ' + str(domain_min[2]) + ' \\\n'
    strings += '    -domain_length             ' + str(domain_length[0]) + ' ' + str(domain_length[1]) + ' ' + str(domain_length[2]) + ' \\\n'
    strings += '    -cfr_steps                 ' + str(cfr_steps[0]) + ' ' + str(cfr_steps[1]) + ' ' + str(cfr_steps[2]) + ' \\\n'
    strings += '    -cfr_flags                 ' + str(cfr_flags[0]) + ' ' + str(cfr_flags[1]) + ' ' + str(cfr_flags[2]) + ' \\\n'
    strings += '    -restart_flags_and_step    ' + str(restart_fs[0]) + ' ' + str(restart_fs[1]) + ' \\\n'
#    strings += ' > ' + logdir + '/logfile.txt 2>&1'
    strings += ' > ' + '../logfile.txt 2>&1'


    # fout #
    f = open(fname, 'w')
    f.write(strings)
    f.close()



if __name__ == '__main__':
    create_pbs_script(fname_pbs)
    create_run(fname_run)

    os.chmod(fname_pbs, 0o777)
    os.chmod(fname_run, 0o777)




#
#
#chmod 777 fname

