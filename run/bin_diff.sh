#!/bin/bash

submit_job () {
    cd $1
    sbatch --reservation daytime $HOME/code/citylbm/run/jd
}

wait_job () {
    while [ -n "$(squeue | awk 'NR>1{print}')" ]; do
        printf "."
        sleep 10;
    done
    echo
}

diff_bin () {
    $HOME/code/citylbm/backup-output/binaries/bin_diff_dirs.py \
        $HOME/code/citylbm/io/iofiles/restart \
        $HOME/code/citylbm-origin/io/iofiles/restart \
        -v
}

submit_job $HOME/code/citylbm
submit_job $HOME/code/citylbm-origin

wait_job

diff_bin | tee ../log/diff.log

