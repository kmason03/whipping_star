#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already
# source mrcnnjob.sh /cluster/tufts/wongjiradlab/jmills09/ubdl_gen2/ /cluster/tufts/wongjiradlab/jmills09/maskrcnn_gen2/ 0
# python tools/save_output_objects.py --dataset particle --cfg configs/tuftscluster_config_0.yaml --load_ckpt weights/u_plane.pth --input_file /cluster/tufts/wongjiradlab/larbys/data/mcc9/mcc9_v29e_dl_run3_G1_extbnb_dlana/data/mcc9_v29e_dl_run3_G1_extbnb_dlana/merged_dlana_d9679e9b-3be3-4411-bc25-6e2cea860827.root --output_dir output2_delete/ --num_images 999
WHIP_DIR=$1
WORK_DIR=$2
SLURM_ARRAY_TASK_ID=$3
export OMP_NUM_THREADS=4
echo "Here"
echo $WHIP_DIR
echo $WORK_DIR
echo $SLURM_ARRAY_TASK_ID
ARRAYID=$SLURM_ARRAY_TASK_ID


cd $WHIP_DIR
source thatroot.sh
cd $WORK_DIR
./DLNuMuDisappearenceSens_OneBinSaveTest --mi ${ARRAYID} >> log_unithrow_shapeonly_${ARRAYID}.txt
