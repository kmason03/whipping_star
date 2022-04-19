#!/bin/bash

#SBATCH --job-name=siginject
#SBATCH --output=log-siginject
#SBATCH --partition wongjiradlab,preempt
#SBATCH --time=0-40:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --array=0-5

CONTAINER=/cluster/tufts/wongjiradlab/larbys/larbys-containers/ubdl_depsonly_py3.6.11_u16.04_cu11_pytorch1.7.1.simg
TOP_DIR=/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/
SCRIPT_DIR=/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/gridscripts/
WORK_DIR=/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/build/Appearence/

module load singularity/3.5.3
cvmfs_config probe fermilab.opensciencegrid.org uboone.opensciencegrid.org
singularity exec -B /cluster:/cluster,/cvmfs:/cvmfs ${CONTAINER} bash -c "cd ${SCRIPT_DIR} && source runjob_fakedata.sh ${TOP_DIR} ${WORK_DIR}"
