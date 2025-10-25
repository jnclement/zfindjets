#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Need tag argument (string), evtnum (int)"
    exit 1
fi
rn=0
evt=$2
nfile=1000
#for i in {0..25000}; do
#    mkdir -p /sphenix/tg/tg01/jets/jocl/evt/$i
#    mkdir -p /sphenix/tg/tg01/jets/jocl/err/$i
#    mkdir -p /sphenix/tg/tg01/jets/jocl/out/$i
#    if [ $(($i % 100)) -eq 0 ]; then
#	echo $i
#    fi
#done
if [ $evt -gt 1000 ]; then
    evt=0
fi
bash run_everysim.sh $1 $nfile $rn 1 $evt 0 $3 $4
