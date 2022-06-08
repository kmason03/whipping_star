#!/bin/bash

#SBATCH --job-name=FCFull
#SBATCH --output=log-FCFull
#SBATCH --partition wongjiradlab
#SBATCH --time=0-60:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --array=923,926,929,1085,1164,1242,1317
#SBATCH --exclude=i2cmp001,i2cmp023

CONTAINER=/cluster/tufts/wongjiradlab/larbys/larbys-containers/ubdl_depsonly_py3.6.11_u16.04_cu11_pytorch1.7.1.simg
TOP_DIR=/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/
SCRIPT_DIR=/cluster/tufts/wongjiradlabnu/kmason03/taritreetests/whipping_star/gridscripts/
WORK_DIR=/cluster/tufts/wongjiradlabnu/kmason03/taritreetests/whipping_star/build/Appearence/


let timer="${SLURM_ARRAY_TASK_ID}*5"
sleep timer

module load singularity/3.5.3
cvmfs_config probe fermilab.opensciencegrid.org uboone.opensciencegrid.org
singularity exec -B /cluster:/cluster,/cvmfs:/cvmfs ${CONTAINER} bash -c "cd ${SCRIPT_DIR} && source runjob_FCFull.sh ${TOP_DIR} ${WORK_DIR}"
