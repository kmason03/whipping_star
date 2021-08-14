#!/bin/bash

function EmptyFolder() 
{
   # remove existing file in the folder
   if [ ! -z "$(ls -A $PWD)" ]; then
     echo -e "Start fresh in the dir..\n"
     rm *
   fi
}

#Some book-keeping, as well as XML locations
TAG="Aug14_UnblindedData_v3_FitA"
COV="/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/FinalSignalBoxUnblinging_Jun2021/Covar/April2021_Full_FitA_v2_FluxXSDetG4.SBNcovar.root"
XML="/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/FinalSignalBoxUnblinging_Jun2021/XML/Redo/master_mc_FITA_final_selection_1g1p_1g0p_2g1p_2g0p_June2021_signal_box_v2.xml"

## SAFE DATA
#DATAXML="/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/FinalSignalBoxUnblinging_Jun2021/XML/master_placeholder_data_final_selection_1g1p_1g0p_2g1p_2g0p_June2021_signal_box_v1.xml"
#DATAXML="/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/FinalSignalBoxUnblinging_Jun2021/XML/master_zero_data_final_selection_1g1p_1g0p_2g1p_2g0p_June2021_signal_box_v1.xml"

## REAL data will be here!
DATAXML="/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/FinalSignalBoxUnblinging_Jun2021/XML/Redo/master_real_FITA_unblinded_data_final_selection_1g1p_1g0p_2g1p_2g0p_June2021_signal_box_v2.xml"

echo -e "\nStarting to run Fit A over real data!!\n"


# directory to copy executables from, and local exetuable name
copy_dir="/uboone/app/users/gge/singlephoton/whipping_star/build/bin/"
fit_exe="sbnfit_singlephoton_unblinding"
constrain_exe="sbnfit_conditional_constraint_unblinding"


## Starting to make necessary obeject
mkdir -p "FitA"
cd "FitA"
EmptyFolder

cp $copy_dir"sbnfit_singlephoton" $fit_exe
cp $copy_dir"sbnfit_conditional_constraint" $constrain_exe

## Run over in conditional constraint

echo -e "First, make Data and MC spectra\n"
./$fit_exe -x $DATAXML -t $TAG"_Data" -m "gen" &> gen_data.log
./$fit_exe -x $XML -t $TAG"_MC" -m "gen" &> gen_mc.log

echo -e "Second, fit!\n"
./$fit_exe -x $XML -t $TAG"_MC" -m "fit" -d $TAG"_Data_CV.SBNspec.root" -c $COV &> fit.log
./$fit_exe -x $XML -t $TAG"_MC" -m "plot" &> plot.log
cp $TAG"_MC_h_delta_chi.pdf" "OUTPUT_RealData_delta_chi.pdf"

LineNumber=$(awk '/Best/{print NR; exit}' fit.log)
((LineNumber+=1))   # the line where best-fit is printed out
LineContent=$(awk "NR==$LineNumber" fit.log)
LineWords=($(echo $LineContent | tr ": " "\n"))
LineLength=${#LineWords[@]}
FittingSuchannel=${LineWords[(($LineLength-2))]}
BFvalue=${LineWords[(($LineLength-1))]}
echo "Third, grab fitting subchannel and BF values: "$FittingSuchannel" "$BFvalue

echo -e "Now, generate constrain plots at BF point\n"
./$constrain_exe -x $XML -t "DataOverlaid" -s $TAG"_MC_CV.SBNspec.root" -d $TAG"_Data_CV.SBNspec.root" -n 2 -c $COV -o -b $FittingSuchannel,$BFvalue &> constrain.log
    
#pdfjam
echo -e "Result saved as:\n\tOUTPUT_RealData_delta_chi.pdf\n\tOUTPUT_DataOverlaid_BF_combined_conditional_constraint.pdf\n"
/uboone/app/users/markrl/useful_scripts/pdfjam/pdfjam/bin/pdfjam DataOverlaid_BF_1g1p_* DataOverlaid_BF_1g0p_* DataOverlaid_Wchi_BF_1g1p_* DataOverlaid_Wchi_BF_1g0p_* --landscape --nup 2x1 --outfile OUTPUT_DataOverlaid_BF_combined_conditional_constraint.pdf
    #/uboone/app/users/markrl/useful_scripts/pdfjam/pdfjam/bin/pdfjam DataOverlaid_BF_1g1p_* DataOverlaid_BF_1g0p_* --landscape --nup 2x1 --outfile OUTPUT_DataOverlaid_BF_combined_conditional_constraint.pdf &> /dev/null
    #/uboone/app/users/markrl/useful_scripts/pdfjam/pdfjam/bin/pdfjam DataOverlaid_Wchi_BF_1g1p_* DataOverlaid_Wchi_BF_1g0p_* --landscape --nup 2x1 --outfile OUTPUT_DataOverlaid_Wchi_BF_combined_conditional_constraint.pdf &> /dev/null

