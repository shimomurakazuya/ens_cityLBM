#!/bin/bash

if [ $# -ne 2 ]; then
    echo "[usage] $0 FROM TO"
    exit
fi

from=$1
to=$2

mkdir -p $to
mv $from/iofiles $from/output $to/
