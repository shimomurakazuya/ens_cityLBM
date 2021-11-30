#!/bin/sh
#$ -cwd
#$ -l f_node=2
#$ -l h_rt=16:00:00
#$ -N oklahoma4m
#$ -e ../log/e
#$ -o ../log/o

rm run.exe
ln -s ../bin/citylbm.bdw_omp run.exe
chmod 755 run.exe

. /etc/profile.d/modules.sh
module purge
module load intel
module load intel-mpi

export OMP_NUM_THREADS=7
export TMPDIR=/tmp

echo $PE_HOSTFILE
cat $PE_HOSTFILE | awk '{}{print $1}' | tee hostfile.txt

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
export citylbm_lxy=4096 # origin 4096 # lxy%32 == 0
export citylbm_nz=$((160/citylbm_dx)) # 160 for dx=1, 80 for dx=2
export citylbm_crft=$((2000/citylbm_dx)) # 60 sec; see --velocity_lbm

ulimit -c unlimited

mpiexec.hydra -f hostfile.txt -genvall -ppn 4 \
       -n $((citylbm_mx * citylbm_my* citylbm_n_ens * citylbm_ngpu_factor)) \
       ./run.exe \
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
|& tee ../log/logfile.txt
