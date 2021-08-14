#!/bin/bash

function EmptyFolder() 
{
   # remove existing file in the folder
   if [ ! -z "$(ls -A $PWD)" ]; then
     echo -e "Start fresh in the dir..\n"
     rm *.log
     rm *.root
   fi
}

## Some book-keeping, as well as XML locations
TAG="Aug14_Unblinding_ExclusiveFit_v3"
COVLoc="/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/Unblinding/WithDuplicate/Exclusive_Fit/gen_covariance_matrix/"

## REAL data will be here!
DATAXML="/uboone/app/users/markrl/SBNfit_uBooNE/July2020_SL7/whipping_star/build/bin/FinalSignalBoxUnblinging_Jun2021/XML/Redo/master_real_unblinded_data_final_selection_1g1p_1g0p_2g1p_2g0p_June2021_signal_box_v2.xml"

## groups of channels to do exclusive fit
ChannelGroup=("1g1p_2g1p" "1g1p_2g1p_2g0p" "1g0p_2g0p" "1g0p_2g1p_2g0p")

## directory to copy executables from, and local exetuable name
copy_dir="/uboone/app/users/gge/singlephoton/whipping_star/build/bin/"
fit_exe="sbnfit_singlephoton_unblinding"
constrain_exe="sbnfit_conditional_constraint_unblinding_poisson"


## OUTdir
OutDir=$PWD


## now, start running

## first, generate data spectra
cp $copy_dir"sbnfit_singlephoton" $fit_exe
echo -e "Generate data spectra\n"
./$fit_exe -x $DATAXML -t $TAG"_Data" -m "gen" &> gen_data.log


# now iterate through different channel groups
for i in ${ChannelGroup[@]};do
    echo -e "Now in directory "$i
    cd $OutDir"/"$i
    #EmptyFolder


    ## copy executables
    cp $copy_dir"sbnfit_singlephoton" $fit_exe
    cp $copy_dir"sbnfit_conditional_constraint" $constrain_exe

    #MC xml
    MCXML="master_mc_final_selection_"$i"_June2021_signal_box_v2.xml"
    COV=$COVLoc"Unblinding_v1_FitE_"$i"_FluxXSDetG4.SBNcovar.root"

    echo -e "start to generate MC spectra\n"
    ./$fit_exe -x $MCXML -t $TAG"_MC" -m "gen" &> gen_mc.log

    echo -e "start to do the fit\n"
    ./$fit_exe -x $MCXML -t $TAG"_MC" -m "fit" -d "../"$TAG"_Data_CV.SBNspec.root" -c $COV &> fit.log
    ./$fit_exe -x $MCXML -t $TAG"_MC" -m "plot" &> plot.log    

    echo -e "Output is saved as OUTPUT_RealData_"$i"_delta_chi.pdf\n"
    cp $TAG"_MC_h_delta_chi.pdf" "OUTPUT_RealData_"$i"_delta_chi.pdf"


    echo -e "Now, start generate conditional constraint CV with data overlaid"
    ./$constrain_exe -x $MCXML -t $TAG"_"$i -s $TAG"_MC_CV.SBNspec.root" -d "../"$TAG"_Data_CV.SBNspec.root" -n 1 -c $COV -o &> constrain.log
done

