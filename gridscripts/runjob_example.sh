#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already

UBDL_DIR=$1
WORKDIR=$2
FILELIST=$3

# change  (+0) to i.e. 2000, if you have already run one set of jobs
let arrayid="$SLURM_ARRAY_TASK_ID+0"
echo $arrayid
#print out file number for log
echo $FILEID
pwd
cd UBDL_DIR
# setup container for ubdl
source setenv_py3.sh
source configure.sh

#get input file name
let line=${arrayid}+1
infile=`sed -n ${line}p ${FILELIST}`
echo $infile

#go to your working directory and make a temp directory for this job
cd ${WORKDIR}
mkdir job_${arrayid}
cd job_${arrayid}
# make a temp copy of vertex file
cp $infile .

echo
ls
pwd
#copy your script to temp directory
# this assumes you already made your cxx code
cp ${WORKDIR}/<cxx object> .

# run your command
./<cxx name> ${infile} .  > ${WORKDIR}/logs/${arrayid}_log.txt

# move outputs to dedicated directory
mv <output file name> ${WORKDIR}/outputs/output_${arrayid}.root
# remove temporary job directory
cd ${WORKDIR}
rm -r job_${arrayid}
echo Finished
