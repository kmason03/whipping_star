#!/bin/bash
#SBATCH --job-name=UniThrow
#SBATCH --output=log_unithrow
#SBATCH --partition batch,preempt
#SBATCH --time=0-24:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --array=0-24
#SBATCH --cpus-per-task=4
# should be 0-24




#
CONTAINER=/cluster/tufts/wongjiradlab/larbys/larbys-containers/ubdl_depsonly_py3.6.11_u16.04_cu11_pytorch1.7.1.simg
WHIP_DIR=/cluster/tufts/wongjiradlab/jmills09/whipping_star
WORK_DIR=/cluster/tufts/wongjiradlab/jmills09/whipping_star/build/NuMuDisappearance

module load singularity/3.5.3
cvmfs_config probe fermilab.opensciencegrid.org uboone.opensciencegrid.org
singularity exec -B /cluster:/cluster,/cvmfs:/cvmfs ${CONTAINER} bash -c "cd ${WHIP_DIR} && source UniThrowjob.sh ${WHIP_DIR} ${WORK_DIR} ${SLURM_ARRAY_TASK_ID}"
