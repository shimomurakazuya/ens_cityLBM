#!/bin/sh
#PJM -L "node=2"
#PJM -L "rscunit=rscunit_ft01"
#PJM -L "rscgrp=dvsmall"
#PJM -L "elapse=10:00"
#PJM --mpi "shape=1"
####PJM --mpi "max-proc-per-node=4"
#PJM --mpi "proc=8"
#PJM -s

rm run.exe
ln -s ../bin/citylbm.fugaku_a64fx_omp run.exe
chmod 755 run.exe

module purge
module load lang/tcsds-1.2.25

export OMP_NUM_THREADS=12
export OMPI_MCA_plm_ple_memory_allocation_policy=bind_local
export TMPDIR=/tmp

mkdir -p ../io/iofiles/restart ../io/output

# Compute problem size
export citylbm_dx=4
export citylbm_n_ens=1 # number of ensemble members;; should be 2*n_ens == f_node
export citylbm_ofs_ens=100 # index offset of ensemble members

# citylbm options
export citylbm_addopt="EXPECT_DX_MESH=$citylbm_dx "
if [[ "$@" =~ "--debug" ]]; then
    export citylbm_addopt+="NO_POSTPROCESS_MONITOR "
    export citylbm_addopt+="NO_IODATA_PARAVIEW "
    #export citylbm_addopt+="NO_IOTHREAD "
    #export citylbm_addopt+="NO_INIT_DECOBOCO "
    #export citylbm_addopt+="NO_IODATA_GZBIN "
fi
export citylbm_addopt+="NO_CITYBOX "
#export citylbm_addopt+="ENSEMBLE_NUDGING "
#export citylbm_addopt+="ENSEMBLE_UVWT "
export citylbm_addopt+="ENSEMBLE_SPIKE "
export citylbm_addopt+="SCALAR_UNIQ "
 
# regression test by binary check. 
# 1: terminate after init
# 2: terminate after init (with restart write/read/write at the end of init)
# 3: t = endtime
# 4: t = endtime (with restart write/read/write at the end of init)
[[ "$@" =~ "--debug" ]] && export citylbm_addopt+="REGRESSION_TEST=3 "

# mpi
export citylbm_ngpu_factor=2 # TSUBAME
export citylbm_mx=$((8/citylbm_dx))
export citylbm_my=$((8/citylbm_dx))
export citylbm_lxy=1024 # origin 4096 # lxy%32 == 0
export citylbm_nz=$((160/citylbm_dx)) # 160 for dx=1, 80 for dx=2
export citylbm_crft=$((40/citylbm_dx)) # 60 sec; see --velocity_lbm

# if n is specified explicitly, it does not work for some reason
#mpiexec -n $((citylbm_mx * citylbm_my* citylbm_n_ens * citylbm_ngpu_factor)) \
mpiexec ./run.exe \
        -n_ensemble_members        $citylbm_n_ens \
        -ofs_ensemble_idx          $citylbm_ofs_ens \
        -gpu_per_node              4 \
        -time_end                  24300.0 \
        -velocity_lbm              10.0  0.100 \
        -number_of_grid_point      $citylbm_nz   2 \
        -domain_min                -$((citylbm_lxy/2))  -$((citylbm_lxy/2))    -8.0  \
        -domain_length              $((citylbm_lxy  ))   $((citylbm_lxy  ))  2560.0  \
        -cfr_steps                 $citylbm_crft $citylbm_crft \
        -cfr_flags                 1  1 \
        -restart_flags_and_step    0 99999 \
        -iofield_freq              10 \
> ../log/logfile.txt 2>&1
