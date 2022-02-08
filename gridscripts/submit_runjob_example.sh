#!/bin/bash

#SBATCH --job-name=example-kmason03
#SBATCH --output=log-example-kmason03
#SBATCH --partition batch
#SBATCH --time=0-2:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-1999

CONTAINER=/cluster/tufts/wongjiradlab/larbys/larbys-containers/ubdl_depsonly_py3.6.11_u16.04_cu11_pytorch1.7.1.simg
UBDL_DIR= <ubdl directory on tufts>
SCRIPT_DIR= <directory of grid scripts>
WORKDIR= <directory code is run>
filelist=/cluster/tufts/wonjiradlab/kmason03/gridscripts/truthdatajobs/bnb_overlay_run3.txt

module load singularity
singularity exec ${CONTAINER} bash -c "cd ${SCRIPT_DIR} && source runjob_example.sh ${UBDL_DIR} ${WORKDIR} ${filelist} "
