#!/bin/bash

TAG=$1
NFILE=$2
RN=$3
ZS=$4
EVT=$5
C2C=$6
if [ $# -lt 6 ]; then
    echo "Need arguments (in order): tag (string), nfile (int), run number (int), ZS (int) EVT (int), C2C (int)"
    exit 1
fi

BASENAME="condor_${TAG}_${NFILE}_${RN}_${ZS}_${EVT}"

SUBNAME="${BASENAME}.sub"

echo "executable = earlysim.sh" > $SUBNAME
echo "arguments = ${TAG} \$(Process) ${RN} ${ZS} ${EVT} ${C2C} ${7}" >> $SUBNAME
echo "output = /sphenix/user/jocl/projects/zfindjets/run/output/out/output_${BASENAME}_\$(Process).out" >> $SUBNAME
echo "request_memory                = 4GB" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "error = /sphenix/user/jocl/projects/zfindjets/run/output/out/output_${BASENAME}_\$(Process).out" >> $SUBNAME
echo "log = /tmp/jocl_${BASENAME}.log" >> $SUBNAME
echo "queue ${NFILE}" >> $SUBNAME

condor_submit $SUBNAME
