#!/bin/bash

#SBATCH --job-name=massgen
#SBATCH --output=log-massgen
#SBATCH --partition batch
#SBATCH --time=0-2:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1

CONTAINER=/cluster/tufts/wongjiradlab/larbys/larbys-containers/ubdl_depsonly_py3.6.11_u16.04_cu11_pytorch1.7.1.simg
TOP_DIR=/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/
SCRIPT_DIR=/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/gridscripts/
WORK_DIR=/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/build/Appearence/

module load singularity/3.5.3
cvmfs_config probe fermilab.opensciencegrid.org uboone.opensciencegrid.org
singularity exec -B /cluster:/cluster,/cvmfs:/cvmfs ${CONTAINER} bash -c "cd ${SCRIPT_DIR} && source runjob_massgen.sh ${TOP_DIR} ${WORK_DIR}"
