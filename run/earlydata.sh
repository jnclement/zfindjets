#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.515
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
export TESTINSTALL=/sphenix/user/jocl/projects/testinstall
echo $LD_LIBRARY_PATH
echo $PATH
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
SUBDIR=${3}
STARTN=$(( $2 * 5 ))
echo $SUBDIR

mkdir -p $SUBDIR
mkdir -p lists
mkdir -p /sphenix/tg/tg01/jets/jocl/zfind/${SUBDIR}/
mkdir -p ./dsts/$SUBDIR
mkdir -p $SUBDIR\_zfind
echo "Made dirs"
cp -r /sphenix/user/jocl/projects/zfindjets/run/run_earlydata.C .
#cp -r /sphenix/user/jocl/projects/zfindjets/run/Fun4All_CaloDataAna.C .
echo "copied root script here"
#cp -r /sphenix/user/jocl/projects/zfindjets/run/lists/dst_jetcalo_run2pp-000$3.list ./lists/$3\_jetcalo.list
echo "copied dstlist here"

cp -r /sphenix/user/jocl/projects/zfindjets/run/lists/dst_calofitting-000$3.list ./lists/$3.list

cp /sphenix/user/jocl/projects/zfindjets/run/lists/dst_triggered_event_seb*-000$3.list ./lists/
echo "copied dstlist here"

for i in {0..4}; do
    UPLN=$(( $STARTN + $i + 1 ))
    #DSTFILE=`sed -n "${UPLN}"p "./lists/${3}_jetcalo.list"`
    #if [ -z "${DSTFILE}" ]; then
    #exit 0
    #fi
    #getinputfiles.pl $DSTFILE
    #FULLPATH=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${DSTFILE}';"`
    #cp $FULLPATH ./$DSTFILE
    #echo "got input files"
    #mv $DSTFILE ./dsts/$3/${3}_$UPLN\_jetcalo.root

    DSTFILE=`sed -n "${UPLN}"p "./lists/${3}.list"`
    if [ -z "${DSTFILE}" ]; then
	exit 0
    fi
    ls dsts
    getinputfiles.pl $DSTFILE
    #FULLPATH=`psql FileCatalog -tA -c "select full_file_path from files where lfn = '${DSTFILE}';"`
    #dd if="${FULLPATH}" of="./${DSTFILE}" bs=12M
    echo "got input files"
    mv $DSTFILE ./dsts/$3/${3}_$UPLN.root
    echo "Running Fun4All now"
    echo $DSTFILE
    echo ./dsts/$3/${3}_$UPLN.root
    root -b -q 'run_earlydata.C("'${1}'",'$UPLN',0,'${5}','${6}','${7}')'
    echo " "
    echo " "
    echo " "
    echo "after run"
    ls -larth
    echo " "
    ls -l $SUBDIR
    ls -l $SUBDIR\_zfind/*
    rm -rf ./dsts/$3/${3}_$UPLN.root
    for file in `ls ./${SUBDIR}_zfind/`; do
	dd if=./${SUBDIR}_zfind/$file of=/sphenix/tg/tg01/jets/jocl/zfind/$SUBDIR/$file
    done
    rm ./${SUBDIR}_zfind/*
done

