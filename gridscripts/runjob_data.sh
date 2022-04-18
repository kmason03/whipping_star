#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already

TOP_DIR=$1
WORK_DIR=$2

# change  (+0) to i.e. 2000, if you have already run one set of jobs
let arrayid="$SLURM_ARRAY_TASK_ID"
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
cp ${WORK_DIR}/DL3plus1_data .

if [ ${arrayid} -eq 0 ]; then
# run your command
./DL3plus1_data 0 >${WORK_DIR}/logs/log_${arrayid}.txt
mv *sinee*txt  ${WORK_DIR}/.
fi

if [ ${arrayid} -eq 1 ]; then
# run your command
./DL3plus1_data 1 >${WORK_DIR}/logs/log_${arrayid}.txt
mv *sinmu*txt  ${WORK_DIR}/.
fi

if [ ${arrayid} -eq 2 ]; then
# run your command
./DL3plus1_data  2 >${WORK_DIR}/logs/log_${arrayid}.txt
mv *sinemu*txt  ${WORK_DIR}/.
fi

# remove temporary job directory
cd ${WORK_DIR}
rm -r job_${arrayid}
echo Finished
