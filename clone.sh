#!/bin/bash -ue

if [ $# -ne 1 ]; then
    echo "Usage: $0 <dst_dir>"
fi

dst=$1

mkdir -p $dst
rsync ./ -auv $dst/

log=$dst/log/clone.log
cat <<EOT > $log
`hostname`
`date`

pwd_origin: $PWD

`git log -n 1`

`git status`

`git diff`
EOT
