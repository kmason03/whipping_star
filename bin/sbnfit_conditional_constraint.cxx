#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TLegend.h"

#include "params.h"
#include "SBNconditional.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

/*************************************************************
 *************************************************************
 *		BEGIN sbnfit_conditional_constraint.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

    std::string xml = "example.xml";

    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"tag", 		required_argument,	0, 't'},
        {"signal", 		required_argument,	0, 's'},
        {"data", 		required_argument,	0, 'd'},
        {"help", 		no_argument,	0, 'h'},
    	{"covar",		required_argument,    0, 'c'},
    	{"genie",		required_argument,    0, 'g'},
        {"flat", required_argument,0,'f'},
        {"zero",no_argument,0,'z'},
        {"cmin",required_argument,0,'k'},
        {"cmax",required_argument,0,'p'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    std::string covar_file = "Stats_Only";
    std::string genie_file = "NONE";
    bool stats_only = true;
    bool use_genie = false;
    std::string signal_file;
    std::string data_file;

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    bool remove_correlations = false;

    double cmin = 0;
    double cmax = -9;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "f:x:s:d:t:c:g:k:p:zh", longopts, &index);

        switch(iarg)
        {
              case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
                break;
              case 'p':
                cmax = (double)strtod(optarg,NULL);
                break;
              case 'k':
                cmin = (double)strtod(optarg,NULL);
                break;
            case 'x':
                xml = optarg;
                break;
            case 'z':
                remove_correlations = true; 
                break;
            case 't':
                tag = optarg;
                break;
     	    case 'c':
    	        covar_file=optarg;
                stats_only = false;
	          break;
	    case 'g':
		use_genie=true;
		genie_file = optarg;
		break;		
            case 's':
                signal_file = optarg;
                break;
	    case 'd':
		data_file = optarg;
		break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_conditional_constraint allows for the plotting of covariance matricies from input root files containing reconstructed variables and covariance matricies. "<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-s\t--signal\t\tInput signal SBNspec.root file"<<std::endl;
                std::cout<<"\t-c\t--covar\t\tInput Systematic Fractional Covariance"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-f\t--flat\t\tAdd a flat percent systematic to fractional covariance matrix (all channels) (default false, pass in percent, i.e 5.0 for 5\% experimental)"<<std::endl;
                std::cout<<"\t-z\t--zero\t\tZero out all off diagonal elements of the systematics covariance matrix (default false, experimental!)"<<std::endl;
                std::cout<<"\t--cmax\t max for fractional covariance plot" << std::endl;
                std::cout<<"\t--cmin\t min for fractional covariance plot" << std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;

                return 0;
        }
    }
    /*************************************************************
     *************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);

    std::cout<<"Begining Covariance Plotting for tag: "<<tag<<std::endl;
    std::cout<<"Loading SBNspec file : "<<signal_file<<" with xml "<<xml<<std::endl;
    SBNspec cv(signal_file,xml, false);
    SBNspec data(data_file,xml, false);
    cv.CollapseVector();
    data.CollapseVector();

    std::vector<SBNspec> vec_spec(2.0, cv);
    std::string constrain_str="NCDeltaRadOverlayLEE";
    //std::string constrain_str="NCDeltaRadOverlaySM";
    //get the NULL spectrum
    vec_spec[0].Scale("NCDeltaRadOverlayLEE", 0.0);
    //vec_spec[0].Scale("NCDeltaRadOverlaySM", 0.0);
    //get the BF values
    vec_spec[1].Scale("NCDeltaRadOverlayLEE", 1.02);
    //vec_spec[1].Scale("NCDeltaRadOverlaySM", 1.92);
    
    std::cout<<"Loading fractional covariance matrix from "<<covar_file<<std::endl;

    TFile * fsys;
    TMatrixT<double> * cov;
    TMatrixT<double>* genie_cov;
    
    if(!stats_only){
        fsys = new TFile(covar_file.c_str(),"read");
        cov = (TMatrixT<double>*)fsys->Get("total_frac_covariance");
	fsys->Close();
   
        TMatrixT<double> frac_flat_matrix(cv.num_bins_total, cv.num_bins_total);
        if(bool_flat_det_sys){
            std::cout << "RUNNING with flat systematics: " << flat_det_sys_percent << "%!" << std::endl;
            frac_flat_matrix.ResizeTo(cv.num_bins_total,cv.num_bins_total);
            frac_flat_matrix.Zero();
            for(int i=0 ; i< cv.num_bins_total; i++){
                   frac_flat_matrix(i,i)=flat_det_sys_percent*flat_det_sys_percent/10000.;
            }
            (*cov) = (*cov)+(frac_flat_matrix);
        }


        if(remove_correlations){
            std::cout<<"WARNING! We are running in   `Remove All Off Diagional Covariances/Correlations Mode` make sure this is what you want. "<<std::endl;
            for(int i=0; i<cv.num_bins_total;i++){ 
                for(int j=0; j<cv.num_bins_total;j++){ 
                    if(i==j)continue;
                    (*cov)(i,j) =0.0;
                }
            }
        }

       int num_nans = 0;
       for(int i=0; i<cv.num_bins_total;i++){ 
                for(int j=0; j<cv.num_bins_total;j++){ 
                    double val = (*cov)(i,j);
                    if( isinf(val) || isnan(val) || val!=val){
                            (*cov)(i,j) = 0.0;
                            num_nans++;
                    }
                }
            }
        std::cout<<"We have removed "<<num_nans<<" / "<<pow(cv.num_bins_total,2)<<" nans in the fractional covariance matrix"<<std::endl;
    }

    //grab genie covariance matrix
    if(use_genie){
	fsys = new TFile(genie_file.c_str(), "read");
	genie_cov  = (TMatrixT<double>*)fsys->Get("individualDir/All_UBGenie_frac_covariance");
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/AxFFCCQEshape_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/DecayAngMEC_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/NormCCCOH_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/NormNCCOH_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/RPA_CCQE_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/Theta_Delta2Npi_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/VecFFCCQEshape_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/XSecShape_CCMEC_UBGenie_frac_covariance"));
	fsys->Close();

	//remove nan's
	for(int i=0; i<cv.num_bins_total;i++){
             for(int j=0; j<cv.num_bins_total;j++){
                    double val = (*genie_cov)(i,j);
                    if( isinf(val) || isnan(val) || val!=val){
                            (*genie_cov)(i,j) = 0.0;
                    }       
             }   
        }  


	*cov -= *genie_cov;

	 //remove correlation in genie
	SBNspec temp_cv(cv);
        temp_cv.Keep(constrain_str, 1.0);
        temp_cv.CalcFullVector();
        std::vector<double> temp_full = temp_cv.full_vector;
        temp_cv=cv;
        temp_cv.Scale(constrain_str, 0.0);
        std::vector<double> temp_others = temp_cv.full_vector;

        for(int i=0; i<temp_full.size(); i++){
           for(int j=0;j<temp_full.size(); j++){
                if(temp_full[i]*temp_full[j]==0 && temp_others[i]*temp_others[j]==0) (*genie_cov)(i,j)=0.0;
           }
        }

	*cov += *genie_cov;
    }

    //starts
    int start_pt = 16;
    int start_pt_1g0p=6;//1g1p only has 6 bins
    

    TFile *fout = new TFile(("Constraint_"+tag+"_output.root").c_str(),"recreate");


    int which_spec_index=0;
    //loop over the spectra vector, and calculating the constrained chi2
    for(SBNspec& fspec:vec_spec){
	    fspec.CalcErrorVector();
	    //if(which_spec_index==0) continue;
	    std::cout << "i = " << which_spec_index<< std::endl;;
	    SBNchi SigChi(fspec, *cov);
	//    SigChi.SetFracPlotBounds(cmin,cmax);
	//    SigChi.PrintMatricies(tag);
	//    sig.WriteOut(tag);


	    double chi_constrain=0;
	    double chi_nueoriginal=0;
	    double chi_total=0;

	    TMatrixT<double> collapsed_covar(cv.num_bins_total_compressed, cv.num_bins_total_compressed);
	    SigChi.FillCollapsedCovarianceMatrix(&collapsed_covar);   //gq: has both nue and numu stats

	    //Add on stats for the constraining bit ONLY
	    for(int i=0; i< cv.num_bins_total_compressed;i++){
		 //add MC intrinsic error
	         collapsed_covar(i,i) += pow(fspec.collapsed_err_vector.at(i), 2.0)- fspec.collapsed_vector.at(i);
	         //if(i >= start_pt) collapsed_covar(i,i) += data.collapsed_vector.at(i);  //use the data stats for constrained matrix
	         if(i >= start_pt) collapsed_covar(i,i) += fspec.collapsed_vector.at(i);  //use the MC stats for constrained matrix
	    }
   
	   double intri_1g1p = 0;
		double overall_1g1p_uncons = 0;
		double intri_1g0p = 0;
		double overall_1g0p_uncons =0;
	   for(int i=0; i<cv.num_bins_total_compressed; i++){
		if(i<start_pt_1g0p) intri_1g1p += pow(fspec.collapsed_err_vector.at(i), 2.0);
		else if(i < start_pt) intri_1g0p+= pow(fspec.collapsed_err_vector.at(i), 2.0);

		for(int j=0 ;j <cv.num_bins_total_compressed; j++){
		  if(i <start_pt_1g0p && j<start_pt_1g0p) overall_1g1p_uncons+=collapsed_covar(i,j);
		  else if( (i>= start_pt_1g0p && i<start_pt) && (j >=start_pt_1g0p && j<start_pt)) overall_1g0p_uncons+=collapsed_covar(i,j);
		}
	  } 
	  


	    // ****************************calculate chi constrain***************************************
	    //everything from here down is collapsed.
	    std::vector<TMatrixT<double>> v_mat =  sbn::splitCovariance(collapsed_covar, start_pt);
	    TMatrixT<double> constrained_mat = sbn::getConstrainedCovariance(v_mat);
	    //constrained nue prediction
	    std::vector<double> constrained_pred = sbn::getConstrainedPrediction(v_mat, fspec.collapsed_vector, data.collapsed_vector, start_pt);

	    TMatrixT<double> full_constrained_mat = constrained_mat;
	    for(int i=0; i<start_pt; i++){
		full_constrained_mat(i,i) +=  ( data.collapsed_vector.at(i) >0.001 ? 3.0/(1.0/data.collapsed_vector.at(i) +  2.0/constrained_pred[i])  : constrained_pred[i]/2.0 );
	    }

	    TMatrixT<double> constrained_invert = SigChi.InvertMatrix(full_constrained_mat);
	    for(int i=0; i< start_pt; i++){
		for(int j=0; j< start_pt; j++){
			chi_constrain += constrained_invert(i,j)*(data.collapsed_vector.at(i)-constrained_pred[i])*(data.collapsed_vector.at(j)-constrained_pred[j]);
		}
	    }
	   
	    double overall_syst_1g1p=0;
	    double overall_syst_1g0p=0;
	    for(int i=0; i< start_pt; i++){
		for(int j=0; j< start_pt; j++){
			if( (i< start_pt_1g0p) && (j<start_pt_1g0p)) overall_syst_1g1p +=constrained_mat(i,j);
			if((i>= start_pt_1g0p) && (j>= start_pt_1g0p)) overall_syst_1g0p +=constrained_mat(i,j);
		}
	    }
 
	    // *************************calculate chi nue original **************************************
	    TMatrixT<double> original_mat = v_mat[0];
	    TMatrixT<double> original_full_mat = original_mat;
	    for(int i=0; i<start_pt; i++){
                original_full_mat(i,i) +=  ( data.collapsed_vector.at(i) >0.001 ? 3.0/(1.0/data.collapsed_vector.at(i) +  2.0/fspec.collapsed_vector.at(i))  : fspec.collapsed_vector.at(i)/2.0 );
            }
	    TMatrixT<double> original_invert = SigChi.InvertMatrix(original_full_mat);
	   for(int i=0; i< start_pt; i++){
                for(int j=0; j< start_pt; j++){
                        chi_nueoriginal += original_invert(i,j)*(data.collapsed_vector.at(i)-fspec.collapsed_vector.at(i))*(data.collapsed_vector.at(j)-fspec.collapsed_vector.at(j));
                }
            }
	    // **************************calculate original nue+numu chi2 ********************************
	    TMatrixT<double> total_collapsed_mat(cv.num_bins_total_compressed, cv.num_bins_total_compressed);
	    total_collapsed_mat = SigChi.FillSystMatrix(*cov, fspec.full_vector, true);
	    TMatrixT<double> total_mat = SigChi.AddStatMatrixCNP(&total_collapsed_mat, fspec.collapsed_vector, data.collapsed_vector);
	    TMatrixT<double> total_invert = SigChi.InvertMatrix(total_mat);
	    for(int i=0; i<cv.num_bins_total_compressed; i++){
		for(int j=0; j< cv.num_bins_total_compressed; j++){
		    chi_total+=total_invert(i,j)*(data.collapsed_vector.at(i)-fspec.collapsed_vector.at(i))*(data.collapsed_vector.at(j)-fspec.collapsed_vector.at(j));
		}
	    }
	    // **************************done calculating chi2 *******************************************
	    std::cout << "Numu-nue side-by-side chi2 value is " << chi_total << std::endl;
	    std::cout << "Nue only original  chi2 value is " << chi_nueoriginal << std::endl;
	    std::cout << "Nue constrained  chi2 value is " << chi_constrain << std::endl;
	    std::cout << "Nue constrained, overall systematic uncertainty is:" << std::endl;
	    std::cout << "\t\t\t\t\t 1g1p: " << sqrt(overall_syst_1g1p) << std::endl;
	    std::cout << "\t\t\t\t\t 1g1p MC intrinsic error " << sqrt(intri_1g1p) << std::endl;
	    std::cout << "\t\t\t\t\t 1g1p overall unconstrained error " << sqrt(overall_1g1p_uncons) << std::endl;
	    std::cout << "\t\t\t\t\t 1g0p MC intrinsic error " << sqrt(intri_1g0p) << std::endl;
	    std::cout << "\t\t\t\t\t 1g0p overall unconstrained error " << sqrt(overall_1g0p_uncons) << std::endl;
	    std::cout << "\t\t\t\t\t 1g0p: " << sqrt(overall_syst_1g0p) << std::endl;
	    std::cout << "constrained 1g1p events: " << std::accumulate(constrained_pred.begin(), constrained_pred.begin()+start_pt_1g0p, 0.0);
	    std::cout << "   constrained 1g1p events: " << std::accumulate(constrained_pred.begin()+start_pt_1g0p, constrained_pred.end(), 0.0)<< std::endl;

	    for(int i=0; i<start_pt; i++){
		std::cout<<i<<" N: "<<fspec.collapsed_vector.at(i)<<" Original: "<<sqrt(collapsed_covar(i,i))/fspec.collapsed_vector.at(i)<<" New: "<<sqrt(constrained_mat(i,i))/fspec.collapsed_vector.at(i)<<" Ratio: "<<sqrt(constrained_mat(i,i))/sqrt(collapsed_covar(i,i))<<std::endl;
	    }



	   // ***************************start drawing histograms******************************************
	   TH1D* h_1g1p_original=(TH1D*)(cv.hist[0]).Clone(); h_1g1p_original->Reset();
	   TH1D* h_1g1p_constrain=(TH1D*)(cv.hist[0]).Clone();  h_1g1p_constrain->Reset();
	   TH1D* h_1g1p_data=(TH1D*)(cv.hist[0]).Clone();  h_1g1p_data->Reset();

	   for(int i=0; i<h_1g1p_data->GetXaxis()->GetNbins(); i++){
		h_1g1p_original->SetBinContent(i+1, fspec.collapsed_vector.at(i));
		h_1g1p_constrain->SetBinContent(i+1, constrained_pred[i]);
		h_1g1p_data->SetBinContent(i+1, data.collapsed_vector.at(i));

		//reset error bars
		h_1g1p_original->SetBinError(i+1, sqrt(original_mat(i,i)));
		h_1g1p_constrain->SetBinError(i+1, sqrt(constrained_mat(i,i)));
		h_1g1p_data->SetBinError(i+1, sqrt(h_1g1p_data->GetBinContent(i+1)));
	   }

	   TH1D* h_1g0p_original=(TH1D*)(cv.hist[10]).Clone(); h_1g0p_original->Reset();
	   TH1D* h_1g0p_constrain=(TH1D*)(cv.hist[10]).Clone();  h_1g0p_constrain->Reset();
	   TH1D* h_1g0p_data=(TH1D*)(cv.hist[10]).Clone();  h_1g0p_data->Reset();

	   for(int i=0; i<h_1g0p_data->GetXaxis()->GetNbins(); i++){
		int local_index = start_pt_1g0p+i;
		h_1g0p_original->SetBinContent(i+1, fspec.collapsed_vector.at(local_index));
		h_1g0p_constrain->SetBinContent(i+1, constrained_pred[local_index]);
		h_1g0p_data->SetBinContent(i+1, data.collapsed_vector.at(local_index));

		//reset error bars
		h_1g0p_original->SetBinError(i+1, sqrt(original_mat(local_index, local_index)));
		h_1g0p_constrain->SetBinError(i+1, sqrt(constrained_mat(local_index, local_index)));
		h_1g0p_data->SetBinError(i+1, sqrt(h_1g0p_data->GetBinContent(i+1)));
	        //std::cout << " " << sqrt(original_mat(local_index, local_index)) << " " << sqrt(full_constrained_mat(local_index, local_index)) << std::endl;
	   }

	  std::vector<TH1D*> vec_hist_original{h_1g1p_original, h_1g0p_original};
	  std::vector<TH1D*> vec_hist_constrain{h_1g1p_constrain, h_1g0p_constrain};
	  std::vector<TH1D*> vec_hist_data{h_1g1p_data, h_1g0p_data};
	  std::vector<std::string> vec_string{"1g1p", "1g0p"};

	  for(int i=0;i<2; i++){
		  TCanvas* c=new TCanvas("c", "c");
		  gStyle->SetOptStat(0);
		  //gStyle->SetErrorX();
		  TLegend* leg=new TLegend(0.7,0.7,0.9,0.9);
		  vec_hist_original[i]->SetLineColor(kMagenta+3);
		  vec_hist_original[i]->SetLineWidth(2);
		  vec_hist_original[i]->SetLineStyle(kDashed);
		  TH1D* h_1g_original_copy= (TH1D*)vec_hist_original[i]->Clone();
		  vec_hist_original[i]->SetFillColorAlpha(kMagenta -10, 0.8);
		  //vec_hist_original[i]->SetFillColor(kMagenta -10);
		  vec_hist_original[i]->SetTitle(Form("%s; Reco shower energy/GeV; Events", vec_string[i].c_str()));
		  //h_1g1p_original->SetFillStyle(4050);  //only useful for TPad
		  vec_hist_original[i]->SetMaximum(std::max(vec_hist_original[i]->GetMaximum(), vec_hist_data[i]->GetMaximum())*1.8);
		  vec_hist_original[i]->SetMinimum(0);

		  vec_hist_constrain[i]->SetLineColor(kGreen+3);
		  vec_hist_constrain[i]->SetLineWidth(2);
		  TH1D* h_1g_constrain_copy = (TH1D*)vec_hist_constrain[i]->Clone();
		  vec_hist_constrain[i]->SetFillColorAlpha(kGreen-7, 0.9);
		  //vec_hist_constrain[i]->SetFillColor(kGreen -10);
		  gStyle->SetHatchesLineWidth(3);
		  vec_hist_constrain[i]->SetTitle(Form("%s; Reco shower energy/GeV; Events", vec_string[i].c_str()));
		  vec_hist_constrain[i]->SetFillStyle(3345);
		  leg->AddEntry(vec_hist_constrain[i], Form("%s Constrained",vec_string[i].c_str()), "LF");

		  //vec_hist_constrain[i]->Draw("E2");
		  vec_hist_original[i]->Draw("E2");
		  h_1g_original_copy->Draw("HIST SAME");
		  vec_hist_constrain[i]->Draw("E2same");
		  h_1g_constrain_copy->Draw("HIST SAME");
		  leg->AddEntry(vec_hist_original[i], Form("%s Orignal", vec_string[i].c_str()), "LF");

		  leg->AddEntry(vec_hist_data[i], Form("%s Toy Data",vec_string[i].c_str()), "ep");

		  leg->Draw();
		  vec_hist_data[i]->SetMarkerStyle(20);
		  vec_hist_data[i]->SetMarkerColor(kBlack);
		  vec_hist_data[i]->SetMarkerSize(1.2);
		  vec_hist_data[i]->SetLineWidth(2);
		  vec_hist_data[i]->SetLineColor(kBlack);
		  vec_hist_data[i]->Draw("E1same");
		  c->Update();
		  if(which_spec_index ==0){
			 //c->SaveAs("CV_1g1p_constrain_comparison.png", "png");
			 c->SaveAs(("CV_"+vec_string[i]+"_constrain_comparison.pdf").c_str(), "pdf");
			 c->Write(Form("CV_%s_comparison",vec_string[i].c_str()));
		  }
		  if(which_spec_index ==1){
			 c->SaveAs(("BF_"+vec_string[i]+"_constrain_comparison.pdf").c_str(), "pdf");
			 c->Write(Form("BF_%s_comparison",vec_string[i].c_str()));
		  }

	  }

	  which_spec_index++;
    }
/*    TFile *fout = new TFile(("Constraint_"+tag+"_output.root").c_str(),"recreate");
    fout->cd();
    constrained_mat.Write("full_covariance");
    fout->Close();

    TH1D h1g1p  = sig.hist[0];  
    TH1D * a1g1p  =(TH1D*)h1g1p.Clone("after");
    TH1D h1g0p  = sig.hist[12];  
    TH1D * a1g0p  =(TH1D*)h1g0p.Clone("after2");

    h1g1p.Clear();
    h1g0p.Clear();

    a1g0p->Clear();
    a1g1p->Clear();


    int i = 0;
    for(i=0; i< h1g1p.GetNbinsX();i++){
            h1g1p.SetBinContent(i+1,0);
            h1g1p.SetBinError(i+1,100*sqrt(collapsed_covar(i,i))/sig.collapsed_vector.at(i));
    }
     for(int j=0; j< h1g0p.GetNbinsX();j++){
            h1g0p.SetBinContent(j+1,0);
            h1g0p.SetBinError(j+1,100*sqrt(collapsed_covar(i,i))/sig.collapsed_vector.at(i));
            i++;
    }
 
    i = 0;
    for(i=0; i< a1g1p->GetNbinsX();i++){
            a1g1p->SetBinContent(i+1,0);
            a1g1p->SetBinError(i+1,100*sqrt(constrained_mat(i,i))/sig.collapsed_vector.at(i));
    }
     for(int j=0; j< a1g0p->GetNbinsX();j++){
            a1g0p->SetBinContent(j+1,0);
            a1g0p->SetBinError(j+1,100*sqrt(constrained_mat(i,i))/sig.collapsed_vector.at(i));
            i++;
    }
    

    TCanvas *c = new TCanvas("c","c",1400,600);
    c->Divide(2,1);
    TPad*p1=(TPad*)c->cd(1);
   
    TLegend l1g1p(0.1,0.8,0.89,0.89);
    l1g1p.SetNColumns(2);
    l1g1p.SetLineColor(kWhite);
    l1g1p.SetLineWidth(0);
    l1g1p.SetFillStyle(0);

    h1g1p.Draw("hist");
    h1g1p.SetFillColor(kRed-7);
    h1g1p.SetLineColor(kRed-7);
    h1g1p.Draw("E2 same");
    h1g1p.SetTitle("1#gamma1p Final Selection");
    h1g1p.GetXaxis()->SetTitle("Reconstructed Shower Energy [GeV]");
    h1g1p.GetYaxis()->SetTitle("Fractional Error [\%]");
    h1g1p.SetFillStyle(3333);
    h1g1p.SetMaximum(60);
    h1g1p.SetMinimum(0);
    
    l1g1p.AddEntry(&h1g1p,"Flux & XS Systematics","F");

    a1g1p->SetFillColor(kRed-7);
    a1g1p->SetLineColor(kRed-7);
    a1g1p->Draw("E2 same");

    l1g1p.AddEntry(a1g1p,"Constrained Systematics","F");
    l1g1p.Draw();
    p1->RedrawAxis();


    TPad*p2 = (TPad*)c->cd(2);
    
    TLegend l1g0p(0.1,0.8,0.89,0.89);
    l1g0p.SetNColumns(2);
    l1g0p.SetLineColor(kWhite);
    l1g0p.SetLineWidth(0);
    l1g0p.SetFillStyle(0);


    h1g0p.Draw("hist");
    h1g0p.SetFillColor(kBlue-7);
    h1g0p.SetLineColor(kBlue-7);
    h1g0p.Draw("E2 same");
    h1g0p.SetTitle("1#gamma0p Final Selection");
    h1g0p.GetXaxis()->SetTitle("Reconstructed Shower Energy [GeV]");
    h1g0p.GetYaxis()->SetTitle("Fractional Error [\%]");
    h1g0p.SetFillStyle(3333);
    h1g0p.SetMaximum(60);
    h1g0p.SetMinimum(0);
    
    a1g0p->SetFillColor(kBlue-7);
    a1g0p->SetLineColor(kBlue-7);
    a1g0p->Draw("E2 same");
    
    l1g0p.AddEntry(&h1g0p,"Flux & XS Systematics","F");
    l1g0p.AddEntry(a1g0p,"Constrained Systematics","F");
    l1g0p.Draw();
    
    p2->RedrawAxis();
    c->SaveAs(("Constraint_"+tag+".pdf").c_str(),"pdf");
*/
    
    fout->Close();
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;
}
