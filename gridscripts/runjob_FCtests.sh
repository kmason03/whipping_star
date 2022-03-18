#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already

TOP_DIR=$1
WORK_DIR=$2

# change  (+0) to i.e. 2000, if you have already run one set of jobs
let arrayid="$SLURM_ARRAY_TASK_ID+0"
echo $arrayid
#print out file number for log
echo $FILEID
pwd
cd ${TOP_DIR}
# setup container for ubdl
source /usr/local/root/root-6.22.06/bin/thisroot.sh

#go to your working directory and make a temp directory for this job
cd ${WORK_DIR}
mkdir job_${arrayid}
cd job_${arrayid}

echo
ls
pwd
#copy your script to temp directory
# this assumes you already made your cxx code
cp ${WORK_DIR}/DL3plus1_FCTests .

# run your command
./DL3plus1_siginject 0  >${WORK_DIR}/logs/log_${arrayid}.txt
ls
mv chis_*txt  ${WORK_DIR}/textfiles/chis_${arrayid}.txt
mv best_*txt  ${WORK_DIR}/textfiles/best_${arrayid}.txt
# remove temporary job directory
cd ${WORK_DIR}
rm -r job_${arrayid}
echo Finished
