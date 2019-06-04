#!/bin/sh

. /etc/profile.d/modules.sh
module load cuda/8.0.61
module load openmpi/2.1.2
module load paraview


qrsh -g jh170031 -l q_core=1 -l h_rt=24:00:00 -pty yes -display $DISPLAY -v TERM /bin/bash
