#!/bin/bash

# sbatch => "Submitted batch job nnn"
get_jid () {
    awk '{}{print $4}'
}

# city
jid_city=`sbatch -p large run/jd-regression --debug | get_jid`
echo "city = $jid_city"

# inflate
cd bindiff/
 jid_inflate=`sbatch -p small -d afterok:$jid_city inflate.sh | get_jid`
 echo "inflate = $jid_inflate"
cd ..

# ens_rankshift
cd bindiff/
 jid_rankshift=`sbatch -p small -d afterany:$jid_inflate ens_rankshift.sh | get_jid`
cd ..

# bindiff
cd bindiff/
 jid_bindiff=`sbatch -p small -d afterany:$jid_rankshift bindiff.sh | get_jid`
 echo "bindiff = $jid_bindiff"
cd ..

tail -F log/logfile.txt bindiff/log_{inflate,ens_rankshift,bindiff}.txt

