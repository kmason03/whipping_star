#!/bin/bash

## script to grab pdf files from gpvms and hack them together

HomeDir="/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/Unblinding"
Script1_path=$HomeDir"/Script1/"
Script3_path=$HomeDir"/Script3/RealData/"



RunOtherScripts=${1:-false}
if [ $RunOtherScripts == "READY" ];then
	scp gge@uboonegpvm04.fnal.gov:/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/FinalSignalBoxUnblinging_Jun2021/LiveUnblinding_24thReal_2Hypo/RES*pdf .
	scp gge@uboonegpvm04.fnal.gov:/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/FinalSignalBoxUnblinging_Jun2021/LiveUnblinding_24thReal_GuanqunCHECK/Script3/RealData/OUTPUT*.pdf .
	scp gge@uboonegpvm04.fnal.gov:/uboone/app/users/markrl/SL7test0/hellstroms_hive/hive/build/bin/Feb2021_EmeraldWorldOrder/1g0p_v2_Unblind/datamc/*RES*pdf .
	scp gge@uboonegpvm04.fnal.gov:/uboone/app/users/markrl/SL7test0/hellstroms_hive/hive/build/bin/Feb2021_EmeraldWorldOrder/1g1p_v2_Unblind/datamc/*RES*pdf .

	pdftk OUTPUT_DataOverlaid_CV_combined_conditional_constraint.pdf RESULT_script2_2Hypo_UNBLINDEDDATA.pdf OUTPUT_RealData_delta_chi.pdf OUTPUT_DataOverlaid_BF_combined_conditional_constraint.pdf RESULT_1g1p_datamc_407series_Jammed.pdf  RESULT_1g0p_datamc_407series_Jammed.pdf cat output OUTPUT_merged.pdf 
else
	rm *.pdf
	scp gge@uboonegpvm04.fnal.gov:$Script1_path"OUTPUT_DataOverlaid_CV_combined_conditional_constraint.pdf" .
fi
