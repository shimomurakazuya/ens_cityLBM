#!/bin/bash -u

#make clean 
make 

if [ "$#" -ne 2 ]; then
    echo "usage: $0 left_dir right_dir"
    exit 1
fi

left=$1/io/iofiles/restart
right=$2/io/iofiles/restart
for d in 0 2 4; do

    files=$( (cd $left/$d; ls {rho,u,v,w,T,scalar}_ptr0_rank*_lv*.dat.gz 2>/dev/null) )
    if [ "$files" != "" ]; then
        ./bindiff.out $left/$d $right/$d $files
    else
        echo restart step $d: not found
    fi

done

