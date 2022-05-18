#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <ctime>
#include <chrono>
// #include <gperftools/profiler.h>

#include "TFile.h"
// #include "TTree.h"
// #include "TH1F.h"
// #include "TString.h"
// #include "TNtuple.h"
// #include "TChain.h"
// #include "TMath.h"
// #include "TSystem.h"
// #include "TMatrixT.h"
// #include "TRandom.h"
// #include "TStyle.h"
// #include "TError.h"
// #include "TCanvas.h"
// #include "TH2F.h"
// #include "TGraph.h"
// #include "TMinuitMinimizer.h"
//#include "TMinuit.h"

// #include "SBNconfig.h"
// #include "SBNchi.h"
// #include "SBNspec.h"
// #include "SBNosc.h"
// #include "SBNfit.h"
// #include "SBNfit3pN.h"
#include "SBNllminimizer.h"
#include "prob.h"
// #include "SBNcovariance.h"


// #define no_argument 0
// #define required_argument 1
// #define optional_argument 2

using namespace sbn;

// define helper functions
// void printbinedges();
// SBNspec GetOscillatedSpectra(const SBNspec& cvSpec, SBNspec massSpec,
// 			float e_app, float e_dis, float m_dis);
TMatrixD GetTotalCov(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Mfracsys);
float GetLLHFromVector(const std::vector<float>& obsSpec, const SBNspec& expSpec,const TMatrixD& Msys,bool prints);
// void testfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
std::string ZeroPadNumber(int num, int digits=3);

// define some global variables
// std::string xml = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/xml/TotalThreePlusOne_full.xml";
// bool draw = true;
// bool printbins = true;
// int mass_start = -1;
// std::string tag = "DL_full";
// // set some start parameters
// const double dm2_lowbound(0.01), dm2_hibound(100);
// const double ue4_lowbound(0.01), ue4_hibound(0.5);
// const double umu4_lowbound(0.01), umu4_hibound(0.5);
// // to genergate dm2_grdpts = 400
// const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
const int nBins_e(12),nBins_mu(19);
const int nBins = nBins_e+nBins_mu;
// const int nFakeExp(1);
// double  mnu, ue, umu;
// int count;

// TMatrixD cov(nBins,nBins);
// SBNspec  appSpec, innerSpec, oscSpec;
// // SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/"+tag+".SBNspec.root",xml);
// SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_Bkg.SBNspec.root",xml);
// std::array<double,nBins>  a_pred;
// std::vector<std::tuple<SBNspec,double>> a_sinsqSpec;
// // TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/build/Appearence/DL_full.SBNcovar_bigbin.root","read");
// // TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins.SBNcovar.root","read");
// TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");

// TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");

int main(int argc, char* argv[]){

  // ProfilerStart("DL3plus1_FCwregen.prof");

  // get the input integer: require 1
  int specific_entry = atoi(argv[1]);
  // int specific_entry = 1;

  // --------------------------------------------------
  // other parameters
  std::string sbnfithome = "/cluster/tufts/wongjiradlabnu/kmason03/taritreetests/whipping_star";
  std::string xml = sbnfithome+"/xml/TotalThreePlusOne_full.xml";
  const int nBins_e(12),nBins_mu(19);
  const int nBins = nBins_e+nBins_mu;
  const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
  const double dm2_lowbound(0.01), dm2_hibound(100);
  const double ue4_lowbound(0.01), ue4_hibound(0.5);
  const double umu4_lowbound(0.01), umu4_hibound(0.5);
  bool draw = true;
  bool printbins = true;
  int mass_start = -1;
  std::string tag = "DL_full";
  // set some start parameters
  // to genergate dm2_grdpts = 400
  const int nFakeExp(1);
  double  mnu, ue, umu;
  int count;

  // --------------------------------------------------

  SBNllminimizer minimizer( xml );

  // make array to  get the parameters
  std::vector<std::vector<float>> params_v;
  // first get the list of points we are throwing universes around
  for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
    for(int uei_base = 0; uei_base < ue4_grdpts; uei_base++){
      for(int umui_base = 0; umui_base < umu4_grdpts; umui_base++){
	std::vector<float> params {mi_base, uei_base, umui_base};
	params_v.push_back(params);
      }
    }
  }

  // open output files
  std::ofstream coordfile;
  std::ofstream chifile;
  std::ofstream covfile;
  std::ofstream spacefile;
  std::ofstream specfile;
  std::string textid = ZeroPadNumber(specific_entry, 5);

  chifile.open("chis_"+textid+".txt", std::ios_base::app);
  covfile.open("cov_tot_big_220314.txt", std::ios_base::app);
  spacefile.open("space.txt", std::ios_base::app);
  specfile.open("spec_FC.txt", std::ios_base::app);
  TFile *fout=new TFile("DLFCTests_hists.root","RECREATE");

  SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_Bkg.SBNspec.root",xml);
  // Load up the necesary bits and store them in vectors on the stack **
  cvSpec.Scale("fullosc",0.0);
  // cvSpec.Scale("1e1p_nue",0.0);
  cvSpec.RemoveMCError();
  std::vector<double> cv_v = cvSpec.collapsed_vector;
  TFile * fsys = new TFile( "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");
  TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");


  SBNchi cvChi(cvSpec, *covFracSys);
  // save cv covar for plotting purposes
  TMatrixD * cvFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
  int a = cvChi.FillCollapsedFractionalMatrix(cvFracSys_collapsed);
  for(short i = 0; i < nBins; i++){
    for(short j = 0; j < nBins; j++){
      covfile<< (*cvFracSys_collapsed)(i,j)<<" ";
    }
    covfile<<std::endl;
  }

  // get the parameters we are testing
  // specific_entry=0;
  int mi_base = params_v[specific_entry][0];
  int uei_base = params_v[specific_entry][1];
  int umui_base = params_v[specific_entry][2];

  // there are 400 mass spectra pregenerated
  int mi_base_new = mi_base*(400/dm2_grdpts);
  //mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
  float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
  float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
  float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
  // calculate scaling factors
  float e_app = 4*pow(ue_base,2)*pow(um_base,2);
  float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
  float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));
  // current test model
  std::cout << "NU Base Model: m41^2:" <<pow(mnu_base,2) <<" m41:"<<mnu_base<< " ue:" << ue_base << " um:" << um_base << std::endl;
  std::cout<<"scale factors: "<<e_app<<" "<<e_dis<<" "<<m_dis<<std::endl;
  NeutrinoModel start_model( mnu_base, uei_base, umui_base );
  start_model.Printall();
  std::cout << "mass tag: " << start_model.mass_tag << std::endl;

  minimizer._gen.regenerate_osc( start_model );
  std::cout << "Full Vector of pre-scale CV spectrum ===========" << std::endl;
  minimizer._gen.spec_central_value.PrintFullVector(true);
  std::cout << "=================================================" << std::endl;
  std::cout << "Full Vector of pre-scale dm2 spectrum ===========" << std::endl;
  minimizer._gen.spec_osc_sinsq.PrintFullVector(true);
  std::cout << "=================================================" << std::endl;

  // get the oscillated spectra
  SBNspec oscSpec =  minimizer.getOscSpectra( mnu_base, ue_base, um_base );
  oscSpec.RemoveMCError();
  //std::vector<double> osc_v = oscSpec.collapsed_vector;
  // for(int i=0;i<nBins;i++){
  // 	specfile<<osc_v[i]<<" ";
  // }
  // specfile<<std::endl;

  if(draw){
    fout->cd();
    TH1D * oscSpec_1e1p_h = new TH1D("oscspec_1e1p","oscspec_1e1p ",nBins_e,0,nBins_e);
    TH1D * oscSpec_1mu1p_h = new TH1D("oscspec_1mu1p","oscspec_1mu1p",nBins_mu,0,nBins_mu);
    std::vector<double> oscspec_v = oscSpec.collapsed_vector;
    for(int i=0;i<nBins;i++){
      if (i<nBins_e) oscSpec_1e1p_h->SetBinContent(i,oscspec_v[i]);
      else oscSpec_1mu1p_h->SetBinContent(i-nBins_e,oscspec_v[i]);
    }
  }

  // I'll be generating a new universe around the oscillated spectrum

  // SBNchi TrueChi(cvSpec, true);
  SBNchi TrueChi(oscSpec, *covFracSys);
  TrueChi.ReloadCoreSpectrum(&oscSpec);
  // int reload =  TrueChi.ReloadCoreSpectrum(*oscSpec);
  fsys->cd();
  TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("collapsed_frac_covariance");
  int b = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);

  // std::cout << "NU Base Model: m41:" <<pow(mnu_base,2) <<" m41^2:"<<mnu_base<< " ue:" << ue_base << " um:" << um_base << std::endl;
  // for(short i = 0; i < nBins; i++){
  //   // for(short j = 0; j < nBins; j++){
  //     std::cout<< (*oscFracSys_collapsed)(i,i)<<std::endl;
  //   // }
  //   // std::cout<<std::endl;
  // }
  std::cout << "CV-SPEC: PrintCollapsedVector =======================" << std::endl;
  cvSpec.PrintCollapsedVector();
  std::cout << "osc-SPEC: PrintFullVector =======================" << std::endl;
  oscSpec.PrintFullVector(true);


  TrueChi.pseudo_from_collapsed = true;
  TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomposition

  // Across several fake experiments for this grid point:
  for(int expi = 0; expi < nFakeExp; expi++){
    std::vector<float> fakeData = TrueChi.GeneratePseudoExperiment();
    // for(int x =0;x<nBins;x++){
    //   std::cout<<fakeData[x]<<" ";
    //   // fakeData[x]=static_cast<int>(fakeData[x]);
    // }
    // std::cout<<std::endl;
    // //
    // fakeData[0] = 3;
    // fakeData[1] = 4;
    // fakeData[2] = 2;
    // fakeData[3] = 7;
    // fakeData[4] = 5;
    // fakeData[5] = 3;
    // fakeData[6] = 3;
    // fakeData[7] = 0;
    // fakeData[8] = 1;
    // fakeData[9] = 2;
    // fakeData[10] = 6;
    // fakeData[11] = 3;
    // fakeData[12] = 30;
    // fakeData[13] = 133;
    // fakeData[14] = 206;
    // fakeData[15] = 254;
    // fakeData[16] = 301;
    // fakeData[17] = 314;
    // fakeData[18] = 399;
    // fakeData[19] = 399;
    // fakeData[20] = 359;
    // fakeData[21] = 363;
    // fakeData[22] = 294;
    // fakeData[23] = 280;
    // fakeData[24] = 208;
    // fakeData[25] = 218;
    // fakeData[26] = 159;
    // fakeData[27] = 93;
    // fakeData[28] = 90;
    // fakeData[29] = 78;
    // fakeData[30] = 51;

    // SBNspec inSpec =  minimizer.getOscSpectra( mnu_base, ue_base, um_base );
		TMatrixD cov_grid =GetTotalCov(fakeData, oscSpec, *oscFracSys_collapsed);
    std::cout<<"here"<<std::endl;
    float gridfit =GetLLHFromVector(fakeData, oscSpec, cov_grid, false);

  	// // start of comptest
  	// // 1 get chi2pt
  	float e_app_in = 4*pow(ue_base,2)*pow(um_base,2);
  	float e_dis_in = 4*pow(ue_base,2)*(1-pow(ue_base,2));
  	float m_dis_in = 4*pow(um_base,2)*(1-pow(um_base,2));


  	//2 get classic grid search
  	float chi_min_grid =1000000;
  	float m41_grid,ue_grid,umu4_grid;

  	// for(int mi_in = 0; mi_in <dm2_grdpts; mi_in++){
  	// 	std::cout<<mi_in<<std::endl;
  	// 	for(int uei_in = 0; uei_in < ue4_grdpts; uei_in++){
  	// 		for(int umui_in = 0; umui_in < umu4_grdpts; umui_in++){
  	// 			float ue_val = pow(10.,(uei_in/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
  	// 			float um_val = pow(10.,(umui_in/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
  	// 			float mnu_val = pow(10.,((mi_in+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
  	// 			// get the oscillated spectra
    //
  	// 			SBNspec inSpec =  minimizer.getOscSpectra( mnu_val, ue_val, um_val );
    //       inSpec.RemoveMCError();
  	// 			SBNchi innerChi(inSpec, *covFracSys);
		// 			TMatrixD * inFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
		// 			int c = innerChi.FillCollapsedFractionalMatrix(inFracSys_collapsed);
  	// 			TMatrixD cov_grid = GetTotalCov(fakeData, inSpec, *inFracSys_collapsed);
    //       // TMatrixD cov_grid = SBNllminimizer::GetTotalCov(fakeData, inSpec, inv_frac_cov );
    //
    //       float chi_test = GetLLHFromVector(fakeData, inSpec, cov_grid, false);
  	// 			delete inFracSys_collapsed;
    //
  	// 			if (chi_test<chi_min_grid){
  	// 				chi_min_grid = chi_test;
  	// 				m41_grid = mnu_val;
  	// 				ue_grid = ue_val;
  	// 				umu4_grid = um_val;
  	// 			}
  	// 		}
  	// 	}
  	// }//end of trad grid inner loops




    double bestfit = minimizer.doFit( fakeData, mnu_base, ue_base, um_base );
    // std::cout<<"here"<<chi_min_grid<<" "<<gridfit<<std::endl;

    // if(bestfit>gridfit) bestfit=gridfit;

    // // save outputs to file
    // chifile<<chi_min_grid<<std::endl;
    chifile<<bestfit<<std::endl;
    chifile<<gridfit<<std::endl;
    //
    // SBNspec testSpec =  minimizer.getOscSpectra( 2.11*2.11, 0.03, 0.24 );
    // SBNchi testChi(testSpec, *covFracSys);
    // TMatrixD * testFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
    // int d = testChi.FillCollapsedFractionalMatrix(testFracSys_collapsed);
    // TMatrixD cov_grid_test =GetTotalCov(fakeData, testSpec, *testFracSys_collapsed);
    // // std::cout<<"test cov: ";
    // // for(short i = 0; i < nBins; i++){
    // //   for(short j = 0; j < nBins; j++){
    // //     std::cout<< cov_grid_test[i][j]<<" ";
    // //   }
    // //   std::cout<<std::endl;
    // // }
    // // std::cout<<"here"<<std::endl;
    // float testfit =GetLLHFromVector(fakeData, testSpec, cov_grid_test, false);
    // std::cout<<"bf pt: "<<testfit<<std::endl;


  } //end of fake experiment loop
// 	std::cout<<mnu_base<<" "<<ue_base<<" "<<um_base<<std::endl;
// 	delete oscFracSys_collapsed;
// // chifile<< std::endl;
// std::cout<<oscFracSys_collapsed->GetNcols()<<std::endl;

  fout->Write();
  fout->Close();

  // ProfilerStop();

  return 0;
} // end of main function

TMatrixD GetTotalCov(const std::vector<float>& obsSpec, const SBNspec& expSpec, const TMatrixD& Mfracsys){
	// function to take the fractional Msys and return totl Msys+Mstat
	// inputs:
	// obsSpec: "data" spectra
	// predSpec: "MC" spectra
	// Mfracsys: fractional (flux+xsec+detvar) covariance matrix
	TMatrixD fullcov(nBins,nBins);
	for(int i = 0; i < nBins; i++){
		for(int j = 0; j < nBins; j++){
			// first set to zero
			fullcov[i][j] = 0.0;
			// scale to the prediction
			fullcov[i][j] = (Mfracsys)[i][j]*expSpec.collapsed_vector[i]*expSpec.collapsed_vector[j];
			// add in stat errors start with CNP for "data" errors
      // std::cout<<fullcov[i][j]<<" "<<obsSpec[i]<<" "<<expSpec.collapsed_vector[i]<<std::endl;
			if(i==j){
				if (expSpec.collapsed_vector[i] >0){
					// 1/observed 2/expectation
					fullcov[i][j] += 3.0 / (1.0/obsSpec[i] + 2.0/expSpec.collapsed_vector[i]);
				}
				else {
					fullcov[i][j] += expSpec.collapsed_vector[i]/2.0;
				}
			}
		}//end of first bin loop
	}//end of second bin loop
	return fullcov;
}//end of GetTotalCov


float GetLLHFromVector(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Msys, bool prints){
	// // function to calculate a chi2 (shape + rate)
	// // inputs:
	// // obsSpec: "data" vector
	// // predSpec: "MC" spectra
	// // Mfracsys: total (flux+xsec+detvar) covariance matrix
  float chisqTest = 0;
  // int nBins = _active_copy->nBins;

  // inv cov for chi2calc
  TMatrixD invcov = Msys;
  invcov.Invert();

  // add the chi2-like part
  chisqTest = 0;
  for(int i = 0; i < nBins; i++){
    for(int j = 0; j < nBins; j++){
// (obsi-predi)*(invcov)*(obsj-predj)
    if(i==j && false) std::cout<<i<<" "
              <<obsSpec[i]<<" "
              <<expSpec.collapsed_vector[i]<<" "
              <<Msys[i][j]<<" "
              <<((obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]))
              <<std::endl;
    chisqTest += (obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]);
        }
      }
  // now need ln(det(2Pi*M))
  // TMatrixD tempcov = 2*3.14159265358979323846*Msys;
  std::cout<<"chi2: "<<chisqTest<<" det: "<<log(Msys.Determinant())<<std::endl;
  chisqTest += log(Msys.Determinant());

  return chisqTest;
}//end of GetLLHFromSpectra

std::string ZeroPadNumber(int num, int digits)
{
  std::stringstream ss;

  // the number is converted to string with the help of stringstream
  ss << num;
  std::string ret;
  ss >> ret;

  // Append zero chars
  int str_length = ret.length();
  for (int i = 0; i < digits - str_length; i++)
    ret = "0" + ret;
  return ret;
}
