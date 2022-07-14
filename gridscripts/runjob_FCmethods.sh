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
cp ${WORK_DIR}/DL3plus1_FCwregen .


# run your command
./DL3plus1_FCwregen 0  > ${WORK_DIR}/logs/${arrayid}_log.txt
mv chis*txt ${WORK_DIR}/tests/pt1/chis_${arrayid}.txt

./DL3plus1_FCwregen 8542  > ${WORK_DIR}/logs/${arrayid}_log.txt
mv chis*txt ${WORK_DIR}/tests/pt2/chis_${arrayid}.txt

./DL3plus1_FCwregen 15624  > ${WORK_DIR}/logs/${arrayid}_log.txt
mv chis*txt ${WORK_DIR}/tests/pt3/chis_${arrayid}.txt

./DL3plus1_FCwregen 9799  > ${WORK_DIR}/logs/${arrayid}_log.txt
mv chis*txt ${WORK_DIR}/tests/pt4/chis_${arrayid}.txt

./DL3plus1_FCwregen 12992  > ${WORK_DIR}/logs/${arrayid}_log.txt
mv chis*txt ${WORK_DIR}/tests/pt5/chis_${arrayid}.txt

./DL3plus1_FCwregen 7370  > ${WORK_DIR}/logs/${arrayid}_log.txt
mv chis*txt ${WORK_DIR}/tests/pt6/chis_${arrayid}.txt

./DL3plus1_FCwregen 9071  > ${WORK_DIR}/logs/${arrayid}_log.txt
mv chis*txt ${WORK_DIR}/tests/pt7/chis_${arrayid}.txt

./DL3plus1_FCwregen 10565  > ${WORK_DIR}/logs/${arrayid}_log.txt
mv chis*txt ${WORK_DIR}/tests/pt8/chis_${arrayid}.txt

# move outputs to dedicated directory
#rm ${WORK_DIR}/textfiles2/chis_${arrayid}.txt 
#mv chis*txt  ${WORK_DIR}/FCresults/.
# remove temporary job directory
cd ${WORK_DIR}
rm -r job_${arrayid}
echo Finished
