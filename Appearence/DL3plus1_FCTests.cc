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
#include "TStyle.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
// #include "TMinuitMinimizer.h"
#include "TMinuit.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"
#include "SBNcovariance.h"
#include "prob.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

// define helper functions
void generate_spectra(bool makecovar, int minrange, int maxrange);
void printbinedges();
SBNspec GetOscillatedSpectra(SBNspec cvSpec, SBNspec massSpec,
			float e_app, float e_dis, float m_dis);
TMatrixD GetTotalCov(std::vector<float>  testSpec, SBNspec predSpec, TMatrixD Mfracsys);
float GetChiSqFromSpectra(std::vector<float>  testSpec, SBNspec predSpec, TMatrixD Msys);
// Int_t testfcn(std::vector<float> obsSpec,TMatrixD Msys,
// 				float m41, float ue4, float um4, bool debug);
void testfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
//float getoscweight(float m41, float ue4, float um4, float ev_e, float ev_L,int type);

// define some global variables
std::string xml = "/uboone/app/users/kmason/whipping_star/xml/TotalThreePlusOne_full.xml";
bool gen = false;
bool printbins = true;
int mass_start = -1;
std::string tag = "DL_full";
// set these parameters at the very start
const double dm2_lowbound(0.01), dm2_hibound(100);
const double ue4_lowbound(0.07), ue4_hibound(0.5);
const double umu4_lowbound(0.07), umu4_hibound(0.5);
// to genergate dm2_grdpts = 100
const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
const int nBins_e(10),nBins_mu(19);
const int nBins = nBins_e+nBins_mu;
const int nFakeExp(1000);
double  mnu, ue, umu;
int count;
std::vector<float> fakeData;
TMatrixD cov(nBins,nBins);
SBNspec  appSpec, innerSpec, oscSpec;
SBNspec cvSpec("MassSpectraData/Full/"+tag+"_Bkg.SBNspec.root",xml);
std::array<double,nBins>  a_pred;
std::ofstream coordfile;
std::ofstream chifile;
std::ofstream covfile;
std::ofstream covtotalfile;
std::ofstream covinvfile;
std::ofstream cvspecfile;
std::ofstream oscspecfile;
Float_t z[5],x[5],y[5],errorz[5];
std::vector<std::tuple<SBNspec,double>> a_sinsqSpec;
TMatrixD * covFracSys_collapsed;

int main(int argc, char* argv[]){
	// start time
	int index;
	int iarg = 0;
	opterr=1;
	auto time_a = std::chrono::steady_clock::now();
	const struct option longopts[] ={
	{"xml", 		required_argument, 	0, 'x'},
	{"gen",	no_argument, 0, 'g'},
	{"part", required_argument,0,'p'},
	{0,			no_argument, 		0,  0},
	};
	// resize fakeData vector
	fakeData.resize(nBins);

	while(iarg != -1){
		iarg = getopt_long(argc,argv, "x:dscp:g", longopts, &index);
		switch(iarg){
			case 'x':
				xml = optarg;
				break;
			case 'g':
				gen = true;
				break;
			case 'p':
				mass_start = atoi(optarg);
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
			return 0;
		}
	}


	//PART 1: precompute all of our sin and sin2 amplitudes so we don't need to later
	if(gen){
		generate_spectra(true, 0,1);
		return 1;
	}

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	if(!gen){

		// open output text files
		coordfile.open("bins_full.txt", std::ios_base::app);
		chifile.open("chis_full_comptest1_full.txt", std::ios_base::app);
		covfile.open("cov_full.txt", std::ios_base::app);
		covtotalfile.open("covtotal_full.txt", std::ios_base::app);
		covinvfile.open("covinv_full.txt", std::ios_base::app);
		cvspecfile.open("cvspec_full.txt", std::ios_base::app);
		oscspecfile.open("oscspec_full.txt", std::ios_base::app);

		// Print binedges for easier plotting
		if(printbins) printbinedges();

		// Load up the necesary bits and store them in vectors on the stack **
		cvSpec.Scale("fullosc",0.0);
		// cvSpec.Scale("extbnb",0.0);
		cvSpec.CollapseVector();
		// cvSpec.PrintFullVector();
		// make a tuple of the spectra and mass term
		for(int mi = 0; mi < 100; mi++){
			mnu = pow(10.,((mi+.5)/float(100)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			std::stringstream stream;
			stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
			std::cout<<mi<<std::endl;
			std::string infile = "MassSpectraData/Full/"+tag+"_SINSQ_dm_"+stream.str()+".SBNspec.root";
			auto inspec = SBNspec(infile,xml);
			// inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
			inspec.CollapseVector();
			std::tuple<SBNspec,double> singletup (inspec,mnu);
			a_sinsqSpec.push_back(singletup);
		}

		// also create the null model for use in minimizer

		//Bring in our covariance matrix!
		// Stats only
		// SBNchi uboone_chi_statsonly(cvSpec,true);

		// Stats + sys
		// Load up cov matrix and add in detector variation component	**
		TFile * fsys = new TFile("katieversion_total.SBNcovar.root","read");
		// TFile * fsys = new TFile("h1_v48_total.SBNcovar.root","read");
		// TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
		// switch to collapsed version
		covFracSys_collapsed = (TMatrixD*)fsys->Get("frac_covariance_collapsed");

		// save covar for plotting purposes
		// for(short i = 0; i < nBins; i++){
		// 	for(short j = 0; j < nBins; j++){
		// 			covfile<< (*covFracSys_collapsed)(i,j)<<" ";
		// 	}
		// 	covfile<<std::endl;
		// }

		cvSpec.CollapseVector();
		cvSpec.PrintCollapsedVector();
		for(int i = 0; i < nBins; i++){
			cvspecfile<<cvSpec.collapsed_vector[i]<<" ";
		}
		cvspecfile<<"\n";

		// first get the central grid point we are throwing universes around
		// for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
		// 	for(int uei_base = 0; uei_base < ue4_grdpts; uei_base++){
		// 		for(int umui_base = 0; umui_base < umu4_grdpts; umui_base++){

		// single, central grid point for testing minimizer
		for(int mi_base = 10; mi_base < 11; mi_base++){
			for(int uei_base = 10; uei_base < 11; uei_base++){
				for(int umui_base = 10; umui_base < 11; umui_base++){
					// there are 100 mass spectra pregenerated
					int mi_base_new = mi_base*(100/dm2_grdpts);
						//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
					float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
					float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
					float mnu_base = pow(10.,((mi_base_new+.5)/float(100)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
					// calculate scaling factors
					float e_app = 4*pow(ue_base,2)*pow(um_base,2);
					float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
					float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));
					// current test model
					std::cout << "NU Base Model: m41^2:" << mnu_base << " ue:" << ue_base << " um:" << um_base << std::endl;
					std::cout<<"scale factors: "<<e_app<<" "<<e_dis<<" "<<m_dis<<std::endl;

					// get the oscillated spectra
					oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
					// oscSpec.PrintCollapsedVector();
					// cvSpec.PrintCollapsedVector();

					 // I'll be generating a new universe around the oscillated spectrum
					// SBNchi TrueChi(cvSpec,*covFracSys);
					// stats only for now to test the rest
					SBNchi TrueChi(cvSpec, true);
					//
					// std::cout<<"made chi"<<std::endl;
					TrueChi.ReloadCoreSpectrum(&oscSpec);
					TrueChi.pseudo_from_collapsed = true;
					TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomposition

					// Across several fake experiments for this grid point:
					for(int expi = 0; expi < nFakeExp; expi++){
						fakeData = TrueChi.GeneratePseudoExperiment();

						// start of comptest
						// 1 get chi2pt
						float e_app_in = 4*pow(ue_base,2)*pow(um_base,2);
						float e_dis_in = 4*pow(ue_base,2)*(1-pow(ue_base,2));
						float m_dis_in = 4*pow(um_base,2)*(1-pow(um_base,2));
						SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app_in, e_dis_in, m_dis_in);
						cov = GetTotalCov(fakeData, inSpec, *covFracSys_collapsed);
						float chi_pT = GetChiSqFromSpectra(fakeData, inSpec, cov);
						chifile<<chi_pT<<" "<<mnu_base<<" "<<ue_base<<" "<<um_base<<std::endl;

						//2 get classic grid search
						float chi_min_grid =1000000;
						float m41_grid,ue_grid,umu4_grid;
						for(int mi_in = 0; mi_in <dm2_grdpts; mi_in++){
							for(int uei_in = 0; uei_in < ue4_grdpts; uei_in++){
								for(int umui_in = 0; umui_in < umu4_grdpts; umui_in++){
									int mi_in_new = mi_in*(100/dm2_grdpts);
									float ue_val = pow(10.,(uei_in/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
									float um_val = pow(10.,(umui_in/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
									float mnu_val = pow(10.,((mi_in_new+.5)/float(100)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
									e_app_in = 4*pow(ue_val,2)*pow(um_val,2);
									e_dis_in = 4*pow(ue_val,2)*(1-pow(ue_val,2));
									m_dis_in = 4*pow(um_val,2)*(1-pow(um_val,2));
									// get the oscillated spectra
									SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_in_new)), e_app_in, e_dis_in, m_dis_in);
									cov = GetTotalCov(fakeData, inSpec, *covFracSys_collapsed);
									float chi_test = GetChiSqFromSpectra(fakeData, inSpec, cov);
									if (chi_test<chi_min_grid){
										chi_min_grid = chi_test;
										m41_grid = mnu_val;
										ue_grid = ue_val;
										umu4_grid = um_val;
									}
								}
							}
						}//end of trad grid inner loops
						chifile<<chi_min_grid<<" "<<m41_grid<<" "<<ue_grid<<" "<<umu4_grid<<std::endl;

						// run the coarse grid search for use later
						float chi_min_gridc =1000000;
						float m41_gridc,ue_gridc,umu4_gridc;
						for(int mi_in = 0; mi_in < 10; mi_in++){
							for(int uei_in = 0; uei_in < 10; uei_in++){
								for(int umui_in = 0; umui_in < 10; umui_in++){
									int mi_in_new = mi_in*(100/dm2_grdpts);
									float ue_val = pow(10.,(uei_in/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
									float um_val = pow(10.,(umui_in/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
									float mnu_val = pow(10.,((mi_in_new+.5)/float(100)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
									e_app_in = 4*pow(ue_val,2)*pow(um_val,2);
									e_dis_in = 4*pow(ue_val,2)*(1-pow(ue_val,2));
									m_dis_in = 4*pow(um_val,2)*(1-pow(um_val,2));
									// get the oscillated spectra
									SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_in_new)), e_app_in, e_dis_in, m_dis_in);
									cov = GetTotalCov(fakeData, inSpec, *covFracSys_collapsed);
									float chi_test = GetChiSqFromSpectra(fakeData, inSpec, cov);
									if (chi_test<chi_min_gridc){
										chi_min_gridc = chi_test;
										m41_gridc = mnu_val;
										ue_gridc = ue_val;
										umu4_gridc = um_val;
									}
								}
							}
						}//end of coarse grid inner loops

						//3. just the minimizer from throw point
						//param #,name, start val,step size, minval,maxval,errtype
						// first some common variables
						Int_t ierflg = 0;
						Double_t dm2_step = 0.01;
						TString dm2_name = "m41";
						Double_t ue4_step = 0.005;
						TString ue4_name = "ue4";
						Double_t umu4_step = 0.005;
						TString umu4_name = "umu4";
						Double_t arglist[2];
						Double_t dm2_val,dm2_err,dm_low,dm_high;
						Int_t dm2_num;
						Double_t ue4_val,ue4_err,ue4_low,ue4_high;
						Int_t ue4_num;
						Double_t umu4_val,umu4_err,umu4_low,umu4_high;
						Int_t umu4_num;
						Int_t npari,nparx,istat;
						// Now ready for minimization step 0 = max steps, 1=tolerance
						// give a large maximum number of steps (tends to take ~125 calls)
						arglist[0] = 1000;
						arglist[1] = 1.;

						// next initialize the minimizer
						TMinuit * gMinuit_simple= new TMinuit(3);
						gMinuit_simple->SetFCN(testfcn);
						Double_t dm2_start = mnu_base;
						Double_t ue4_start = ue_base;
						Double_t umu4_start = um_base;
						// start values near lower bound, but some wiggle room
						gMinuit_simple->mnparm(0,dm2_name,dm2_start,dm2_step,dm2_lowbound,dm2_hibound,ierflg);
						gMinuit_simple->mnparm(1,ue4_name,ue4_start,ue4_step,ue4_lowbound,ue4_hibound,ierflg);
						gMinuit_simple->mnparm(2,umu4_name,umu4_start,umu4_step,umu4_lowbound,umu4_hibound,ierflg);
						gMinuit_simple->mnexcm("MIGRAD", arglist ,2,ierflg);
						Double_t fmin,fedm,errdef;
						gMinuit_simple->mnpout(0,dm2_name,dm2_val,dm2_err,dm_low,dm_high,dm2_num);
						gMinuit_simple->mnpout(1,ue4_name,ue4_val,ue4_err,ue4_low,ue4_high,ue4_num);
						gMinuit_simple->mnpout(2,umu4_name,umu4_val,umu4_err,umu4_low,umu4_high,umu4_num);
						gMinuit_simple->mnstat(fmin,fedm,errdef,npari,nparx,istat);
						chifile <<fmin<<" "<<dm2_val<<" "<<ue4_val<<" "<<umu4_val<<std::endl;
						delete gMinuit_simple;

						// 4. coarse grid start
						// next initialize the minimizer
						TMinuit * gMinuit_coarse= new TMinuit(3);
						gMinuit_coarse->SetFCN(testfcn);
						dm2_start = m41_gridc;
						ue4_start = ue_gridc;
						umu4_start = umu4_gridc;
						// start values near lower bound, but some wiggle room
						gMinuit_coarse->mnparm(0,dm2_name,dm2_start,dm2_step,dm2_lowbound,dm2_hibound,ierflg);
						gMinuit_coarse->mnparm(1,ue4_name,ue4_start,ue4_step,ue4_lowbound,ue4_hibound,ierflg);
						gMinuit_coarse->mnparm(2,umu4_name,umu4_start,umu4_step,umu4_lowbound,umu4_hibound,ierflg);
						gMinuit_coarse->mnexcm("MIGRAD", arglist ,2,ierflg);
						Double_t fminc,fedmc,errdefc;
						gMinuit_coarse->mnpout(0,dm2_name,dm2_val,dm2_err,dm_low,dm_high,dm2_num);
						gMinuit_coarse->mnpout(1,ue4_name,ue4_val,ue4_err,ue4_low,ue4_high,ue4_num);
						gMinuit_coarse->mnpout(2,umu4_name,umu4_val,umu4_err,umu4_low,umu4_high,umu4_num);
						gMinuit_coarse->mnstat(fminc,fedmc,errdefc,npari,nparx,istat);
						chifile <<fminc<<" "<<dm2_val<<" "<<ue4_val<<" "<<umu4_val<<std::endl;
						delete gMinuit_coarse;

						// 5. fine grid start
						// next initialize the minimizer
						TMinuit * gMinuit_fine= new TMinuit(3);
						gMinuit_fine->SetFCN(testfcn);
						dm2_start = m41_grid;
						ue4_start = ue_grid;
						umu4_start = umu4_grid;
						// start values near lower bound, but some wiggle room
						gMinuit_fine->mnparm(0,dm2_name,dm2_start,dm2_step,dm2_lowbound,dm2_hibound,ierflg);
						gMinuit_fine->mnparm(1,ue4_name,ue4_start,ue4_step,ue4_lowbound,ue4_hibound,ierflg);
						gMinuit_fine->mnparm(2,umu4_name,umu4_start,umu4_step,umu4_lowbound,umu4_hibound,ierflg);
						gMinuit_fine->mnexcm("MIGRAD", arglist ,2,ierflg);
						Double_t fminf,fedmf,errdeff;
						gMinuit_fine->mnpout(0,dm2_name,dm2_val,dm2_err,dm_low,dm_high,dm2_num);
						gMinuit_fine->mnpout(1,ue4_name,ue4_val,ue4_err,ue4_low,ue4_high,ue4_num);
						gMinuit_fine->mnpout(2,umu4_name,umu4_val,umu4_err,umu4_low,umu4_high,umu4_num);
						gMinuit_fine->mnstat(fminf,fedmf,errdeff,npari,nparx,istat);
						chifile <<fminf<<" "<<dm2_val<<" "<<ue4_val<<" "<<umu4_val<<std::endl;
						delete gMinuit_fine;



						// start of chi2 space drawing
						// // loop over m (keep u's constant)
						// for(int mi_in = 0; mi_in <100; mi_in++){
						// 	float mnu_test = pow(10.,((mi_in+.5)/float(100)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
						// 	float e_app_in = 4*pow(ue_base,2)*pow(um_base,2);
						// 	float e_dis_in = 4*pow(ue_base,2)*(1-pow(ue_base,2));
						// 	float m_dis_in = 4*pow(um_base,2)*(1-pow(um_base,2));
						// 	SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_in)), e_app_in, e_dis_in, m_dis_in);
						// 	cov = GetTotalCov(fakeData, inSpec, *covFracSys_collapsed);
						// 	float chi_test = GetChiSqFromSpectra(fakeData, inSpec, cov);
						// 	chifile<<chi_test<<" ";
						//
						// }// end of loop over m
						// chifile<<std::endl;
						// std::cout<<"finished loop1"<<std::endl;
						// // loop over ue4 (keep others constant)
						// for(int uei_in = 0; uei_in <100; uei_in++){
						// 	float ue_val = pow(10.,(uei_in/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
						// 	float e_app_in = 4*pow(ue_val,2)*pow(um_base,2);
						// 	float e_dis_in = 4*pow(ue_val,2)*(1-pow(ue_val,2));
						// 	float m_dis_in = 4*pow(um_base,2)*(1-pow(um_base,2));
						// 	SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app_in, e_dis_in, m_dis_in);
						// 	cov = GetTotalCov(fakeData, inSpec, *covFracSys_collapsed);
						// 	float chi_test = GetChiSqFromSpectra(fakeData, inSpec, cov);
						// 	chifile<<chi_test<<" ";
						//
						// }// end of loop over m
						// chifile<<std::endl;
						// std::cout<<"finished loop2"<<std::endl;
						// // loop over ue4 (keep others constant)
						// for(int umui_in = 0; umui_in <100; umui_in++){
						// 	float um_val = pow(10.,(umui_in/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
						// 	float e_app_in = 4*pow(ue_base,2)*pow(um_val,2);
						// 	float e_dis_in = 4*pow(ue_base,2)*(1-pow(ue_base,2));
						// 	float m_dis_in = 4*pow(um_val,2)*(1-pow(um_val,2));
						// 	SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app_in, e_dis_in, m_dis_in);
						// 	cov = GetTotalCov(fakeData, inSpec, *covFracSys_collapsed);
						// 	float chi_test = GetChiSqFromSpectra(fakeData, inSpec, cov);
						// 	chifile<<chi_test<<" ";
						//
						// }// end of loop over m
						// chifile<<std::endl;
						// std::cout<<"finished loop3"<<std::endl;


						// // do an inner loop for a traditional grid search to test:
						// float chi_min_grid = 1000000000000;
						// float m41_grid, ue_grid, umu4_grid;
						// for(int mi_in = 0; mi_in <dm2_grdpts; mi_in++){
						// 	for(int uei_in = 0; uei_in < ue4_grdpts; uei_in++){
						// 		for(int umui_in = 0; umui_in < umu4_grdpts; umui_in++){
						// 			int mi_in_new = mi_in*(100/dm2_grdpts);
						// 			float ue_val = pow(10.,(uei_in/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
						// 			float um_val = pow(10.,(umui_in/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
						// 			float mnu_val = pow(10.,((mi_in_new+.5)/float(100)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
						// 			float e_app_in = 4*pow(ue_val,2)*pow(um_val,2);
						// 			float e_dis_in = 4*pow(ue_val,2)*(1-pow(ue_val,2));
						// 			float m_dis_in = 4*pow(um_val,2)*(1-pow(um_val,2));
						// 			// get the oscillated spectra
						// 			SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_in_new)), e_app_in, e_dis_in, m_dis_in);
						// 			cov = GetTotalCov(fakeData, inSpec, *covFracSys);
						// 			float chi_test = GetChiSqFromSpectra(fakeData, inSpec, cov);
						// 			if (chi_test<chi_min_grid){
						// 				chi_min_grid = chi_test;
						// 				m41_grid = mnu_val;
						// 				ue_grid = ue_val;
						// 				umu4_grid = um_val;
						// 			}
						// 		}
						// 	}
						// }//end of trad grid inner loops



						// // start with a few different stepsizes -- test
						// std::vector<double> usteps {0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,.25};
						// std::vector<double> msteps {0.00001,0.0001,0.001,0.005,0.01,0.05,0.1,.25,0.5,0.75,1.0,2.5,5.,10.,15.};
						//
						//
						// for(int mi_step = 0; mi_step < msteps.size(); mi_step++){
						// 	for(int uei_step = 0; uei_step < usteps.size(); uei_step++){
						// 		for(int umui_step = 0; umui_step < usteps.size(); umui_step++){
						//
						// 			float ue_val = ue_base;
						// 			float um_val = um_base;
						// 			float mnu_val = mnu_base;
						//
						// 			// use these as new minimizer starts
						// 			//param #,name, start val,step size, minval,maxval,errtype
						// 			Int_t ierflg = 0;
						// 			Double_t dm2_step = msteps[mi_step];
						// 			Double_t dm2_start = mnu_val;
						// 			TString dm2_name = "m41";
						// 			Double_t ue4_step = usteps[uei_step];
						// 			Double_t ue4_start = ue_val;
						// 			TString ue4_name = "ue4";
						// 			Double_t umu4_step = usteps[umui_step];
						// 			Double_t umu4_start = um_val;
						// 			TString umu4_name = "umu4";
						//
						//
						// 			// next initialize the minimizer
						// 			TMinuit * gMinuit= new TMinuit(3);
						// 			gMinuit->SetFCN(testfcn);
						//
						// 			// start values near lower bound, but some wiggle room
						// 			gMinuit->mnparm(0,dm2_name,dm2_start,dm2_step,dm2_lowbound,dm2_hibound,ierflg);
						// 			gMinuit->mnparm(1,ue4_name,ue4_start,ue4_step,ue4_lowbound,ue4_hibound,ierflg);
						// 			gMinuit->mnparm(2,umu4_name,umu4_start,umu4_step,umu4_lowbound,umu4_hibound,ierflg);
						// 			//
						// 			Double_t arglist[2];
						// 			// Now ready for minimization step 0 = max steps, 1=tolerance
						// 			// give a large maximum number of steps (tends to take ~125 calls)
						// 			arglist[0] = 1000;
						// 			arglist[1] = 1.;
						// 			gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
						//
						// 			// retrieve final parameters and best fit value
						// 			// probable won't need all these output products,
						// 			// but need to initialize to use call function which fills them
						// 			Double_t dm2_val,dm2_err,dm_low,dm_high;
						// 			Int_t dm2_num;
						// 			Double_t ue4_val,ue4_err,ue4_low,ue4_high;
						// 			Int_t ue4_num;
						// 			Double_t umu4_val,umu4_err,umu4_low,umu4_high;
						// 			Int_t umu4_num;
						// 			Double_t fmin,fedm,errdef;
						// 			Int_t npari,nparx,istat;
						// 			gMinuit->mnpout(0,dm2_name,dm2_val,dm2_err,dm_low,dm_high,dm2_num);
						// 			gMinuit->mnpout(1,ue4_name,ue4_val,ue4_err,ue4_low,ue4_high,ue4_num);
						// 			gMinuit->mnpout(2,umu4_name,umu4_val,umu4_err,umu4_low,umu4_high,umu4_num);
						// 			gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
						//
						// 			cov = GetTotalCov(fakeData, oscSpec, *covFracSys_collapsed);
						// 			// save covar for plotting purposes
						// 			for(short i = 0; i < nBins; i++){
						// 				for(short j = 0; j < nBins; j++){
						// 						covfile<< (cov)(i,j)<<" ";
						// 				}
						// 				covfile<<std::endl;
						// 			}
						//
						// 			float chi_pT = GetChiSqFromSpectra(fakeData, oscSpec, cov);
						// 			chifile <<fmin<<" "<<chi_pT<<" "<<dm2_val<<" "<<ue4_val<<" "<<umu4_val<<" ";
						// 			chifile <<msteps[mi_step]<<" "<<usteps[uei_step]<<" "<<usteps[umui_step]<<" "<<expi<<std::endl;
						// 			delete gMinuit;
						// 		}//end of inner mass step loop
						// 	}//end of inner ue4 step loop
						// }//end of inner umu4 step loop
						std::cout<<mnu_base<<" "<<ue_base<<" "<<um_base<<std::endl;
						// std::cout<<"------------------------------"<<std::endl;
						// std::cout<<"grid search minimum: "<<chi_min_grid<<std::endl;
						// std::cout<<m41_grid<<" "<<ue_grid<<" "<<umu4_grid<<std::endl;
						// std::cout<<"------------------------------"<<std::endl;


						// ------------------------------------------------------------

						// //start with a few different seeds test
						// std::vector<int> ids={0,5,10,15,20,24};
						//
						// for(int mi_in = 0; mi_in < ids.size(); mi_in++){
						// 	for(int uei_in = 0; uei_in < ids.size(); uei_in++){
						// 		for(int umui_in = 0; umui_in < ids.size(); umui_in++){
						// 			// get the parameter value based on id
						// 			int mi_idx = ids[mi_in];
						// 			int uei_idx = ids[uei_in];
						// 			int umui_idx = ids[umui_in];
						// 			float ue_val = pow(10.,(uei_idx/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
						// 			float um_val = pow(10.,(umui_idx/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
						// 			float mnu_val = pow(10.,((mi_idx+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
						//
						// 			// use these as new minimizer starts
						// 			//param #,name, start val,step size, minval,maxval,errtype
						// 			Int_t ierflg = 0;
						// 			Double_t dm2_step = 1;
						// 			Double_t dm2_start = mnu_val;
						// 			TString dm2_name = "m41";
						// 			Double_t ue4_step = .00001;
						// 			Double_t ue4_start = ue_val;
						// 			TString ue4_name = "ue4";
						// 			Double_t umu4_step = .00001;
						// 			Double_t umu4_start = um_val;
						// 			TString umu4_name = "umu4";
						//
						// 			// next initialize the minimizer
						// 			TMinuit * gMinuit= new TMinuit(3);
						// 			gMinuit->SetFCN(testfcn);
						//
						// 			// start values near lower bound, but some wiggle room
						// 			gMinuit->mnparm(0,dm2_name,dm2_start,dm2_step,dm2_lowbound,dm2_hibound,ierflg);
						// 			gMinuit->mnparm(1,ue4_name,ue4_start,ue4_step,ue4_lowbound,ue4_hibound,ierflg);
						// 			gMinuit->mnparm(2,umu4_name,umu4_start,umu4_step,umu4_lowbound,umu4_hibound,ierflg);
						// 			//
						// 			Double_t arglist[2];
						// 			// Now ready for minimization step 0 = max steps, 1=tolerance
						// 			// give a large maximum number of steps (tends to take ~125 calls)
						// 			arglist[0] = 1000;
						// 			arglist[1] = 1.;
						// 			gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
						//
						// 			// retrieve final parameters and best fit value
						// 			// probable won't need all these output products,
						// 			// but need to initialize to use call function which fills them
						// 			Double_t dm2_val,dm2_err,dm_low,dm_high;
						// 			Int_t dm2_num;
						// 			Double_t ue4_val,ue4_err,ue4_low,ue4_high;
						// 			Int_t ue4_num;
						// 			Double_t umu4_val,umu4_err,umu4_low,umu4_high;
						// 			Int_t umu4_num;
						// 			Double_t fmin,fedm,errdef;
						// 			Int_t npari,nparx,istat;
						// 			gMinuit->mnpout(0,dm2_name,dm2_val,dm2_err,dm_low,dm_high,dm2_num);
						// 			gMinuit->mnpout(1,ue4_name,ue4_val,ue4_err,ue4_low,ue4_high,ue4_num);
						// 			gMinuit->mnpout(2,umu4_name,umu4_val,umu4_err,umu4_low,umu4_high,umu4_num);
						// 			gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
						//
						// 			// some cout statements to test
						//
						// 			cov = GetTotalCov(fakeData, oscSpec, *covFracSys_collapsed);
						// 			// save covar for plotting purposes
						// 			for(short i = 0; i < nBins; i++){
						// 				for(short j = 0; j < nBins; j++){
						// 						covfile<< (cov)(i,j)<<" ";
						// 				}
						// 				covfile<<std::endl;
						// 			}
						//
						// 			float chi_pT = GetChiSqFromSpectra(fakeData, oscSpec, cov);
						// 			// see what the chi2 is at the coarse grid pt.
						// 			float e_app_in = 4*pow(ue_val,2)*pow(um_val,2);
						// 			float e_dis_in = 4*pow(ue_val,2)*(1-pow(ue_val,2));
						// 			float m_dis_in = 4*pow(um_val,2)*(1-pow(um_val,2));
						// 			// get the oscillated spectra
						// 			SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_in)), e_app_in, e_dis_in, m_dis_in);
						// 			cov = GetTotalCov(fakeData, inSpec, *covFracSys_collapsed);
						// 			float chi_coarse = GetChiSqFromSpectra(fakeData, inSpec, cov);
						//
						//
						// 			// std::cout<<"chi pt: "<<chi_pT<<std::endl;
						// 			chifile <<fmin<<" "<<chi_pT<<" "<<dm2_val<<" "<<ue4_val<<" "<<umu4_val<<" "<<chi_coarse<<std::endl;
						// 			// chifile <<chi_pT-fmin;
						// 			std::cout<<mnu_base<<" "<<ue_base<<" "<<um_base<<std::endl;
						// 			delete gMinuit;
						// 		}//end of inner mass loop
						// 	}//end of inner ue4 loop
						// }//end of inner umu4 loop
						//----------------------------------------------------------------

						// // do an inner loop for a traditional grid search to test:
						// float chi_min_grid = 10000000;
						// float m41_grid, ue_grid, umu4_grid;
						// for(int mi_in = 0; mi_in <dm2_grdpts; mi_in++){
						// 	for(int uei_in = 0; uei_in < ue4_grdpts; uei_in++){
						// 		for(int umui_in = 0; umui_in < umu4_grdpts; umui_in++){
						// 			float ue_val = pow(10.,(uei_in/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
						// 			float um_val = pow(10.,(umui_in/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
						// 			float mnu_val = pow(10.,((mi_in+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
						// 			float e_app_in = 4*pow(ue_val,2)*pow(um_val,2);
						// 			float e_dis_in = 4*pow(ue_val,2)*(1-pow(ue_val,2));
						// 			float m_dis_in = 4*pow(um_val,2)*(1-pow(um_val,2));
						// 			// get the oscillated spectra
						// 			SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_in)), e_app_in, e_dis_in, m_dis_in);
						// 			cov = GetTotalCov(fakeData, inSpec, *covFracSys_collapsed);
						// 			float chi_test = GetChiSqFromSpectra(fakeData, inSpec, cov);
						// 			if (chi_test<chi_min_grid){
						// 				chi_min_grid = chi_test;
						// 				m41_grid = mnu_val;
						// 				ue_grid = ue_val;
						// 				umu4_grid = um_val;
						// 			}
						// 		}
						// 	}
						// }//end of trad grid inner loops
						// std::cout<<"------------------------------"<<std::endl;
						// std::cout<<"grid search minimum: "<<chi_min_grid<<std::endl;
						// std::cout<<m41_grid<<" "<<ue_grid<<" "<<umu4_grid<<std::endl;
						// std::cout<<"------------------------------"<<std::endl;


					} //end of universeloop
					std::cout<<mnu_base<<" "<<ue_base<<" "<<um_base<<std::endl;
				// chifile<< std::endl;
				}// end of loop over base umu
			}//end of loop over base ue
		}//end of loop over base mass
	}//end of not gen loop

	return 0;
} // end of main function


//-----------------------------------------------------------------------------
//----------------------HELPER FUNCTIONS---------------------------------------
//----------------------------------------------------------------------------

 void generate_spectra(bool makecovar, int minrange, int maxrange){
	 //function to generate the different mass spectra we need
	 // prerunning this speeds up the initial grid search part and helps with plotting
	 // note: currently run twice to prevent the job being killed
	 // inputs:
	 // bool makecovar: true if new covar needed
	 // int minrange: min mass index to run
	 // int maxrange: max mass index to run
	 // outputs:
	 // none, but writes spectra in root files to current directory
	if(makecovar){
			SBNcovariance _covar(xml);
		 _covar.FormCovarianceMatrix(tag);
		 _covar.PrintMatricies(tag);
		 _covar.frac_covariance.Print();
	}
	// start with a null model
	NeutrinoModel nullModel(0, 0, 0);
	SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
	SBNspec bkg = bkgo->spec_central_value;
	bkg.WriteOut(tag+"_Bkg");

	float mnu;
	for(int mi = minrange; mi < maxrange; mi++){
		mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		// float stepsize= (dm2_hibound-dm2_lowbound)/float(dm2_grdpts);
		// mnu=dm2_lowbound+mi*stepsize;

		//Model: mnu, ue4, um4
		//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
		//see davio's thesis eqs 4.14-16
		//we're precomputing, so we don't really care about the u's set to sin2=1
		NeutrinoModel testModel(mnu, sqrt(.5), sqrt(.5));

		// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
		SBNgenerate * gen = new SBNgenerate(xml,testModel);

		// Write them to files
		gen->WritePrecomputedOscSpecs(tag);
	}
	return;
}//end of generate spectra function

void printbinedges(){
	// funtion that prints the bins to the output textfile
	for(int mi = 0; mi < 100; mi++){
		mnu = pow(10.,((mi)/float(100)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		coordfile << pow(mnu,2) << " ";
	}
	coordfile << std::endl;

	for(int uei = 0; uei <100; uei++){
		ue = pow(10.,(uei/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
		coordfile << ue << " ";
	}
	coordfile << std::endl;

	for(int umui = 0; umui < 100; umui++){
		umu = pow(10.,(umui/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
		coordfile << umu << " ";
	}
	coordfile << std::endl;
	return;
}//end of print bins function


SBNspec GetOscillatedSpectra(SBNspec cvSpec, SBNspec massSpec,
			float e_app, float e_dis, float m_dis){
	// function to take the cv spec and return the oscillated spec
	// this is only set up to do all osc at once.
	// inputs:
	// cvspec: cv spectum
	// massspec: oscillated spec for the given mass, maximum oscillations
	// e_app,e_dis,m_dis oscillation scaling parameters
	// output oscspec: oscillated spectrum
	SBNspec testSpec = cvSpec;
	testSpec.Scale("fullosc",0.0);
	massSpec.Scale("fullosc",e_app);
	massSpec.Scale("bnb",-1*m_dis);
	massSpec.Scale("nue",-1*e_dis);
	massSpec.Scale("ccpi0",-1*m_dis);
	massSpec.Scale("ncpi0",-1*m_dis);
	massSpec.Scale("ext",0.0);
	// massSpec.PrintCollapsedVector();
	testSpec.Add(&massSpec);
	// std::cout<<"oscillated spec"<<std::endl;
	// testSpec.PrintCollapsedVector();
	return testSpec;
}//end of GetOscillatedSpectra

TMatrixD GetTotalCov(std::vector<float>  testSpec, SBNspec predSpec, TMatrixD Mfracsys){
	// function to take the fractional Msys and return total Msys+Mstat
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
			fullcov[i][j] = (Mfracsys)[i][j]*predSpec.collapsed_vector[i]*predSpec.collapsed_vector[j];
			// add in stat errors start with CNP for "data" errors
			if(i==j){
				if (predSpec.collapsed_vector[i] >0 ){
					fullcov[i][j] += 3.0 / (1.0/testSpec[i] + 2.0/predSpec.collapsed_vector[i]);
					// cov[i][j] += 3.0 / (1.0/testSpec.collapsed_vector[i] + 2.0/predSpec.collapsed_vector[i]);
				}
				else {
					fullcov[i][j] += predSpec.collapsed_vector[i]/2.0;
				}
			}
		}//end of first bin loop
	}//end of second bin loop
	return fullcov;
}//end of GetTotalCov

float GetChiSqFromSpectra(std::vector<float>  testSpec, SBNspec predSpec, TMatrixD Msys){
	// function to calculate a chi2 (shape + rate)
	// inputs:
	// obsSpec: "data" spectra
	// predSpec: "MC" spectra
	// Mfracsys: total (flux+xsec+detvar) covariance matrix
	float chisqTest;

	// inv cov for chi2calc
	TMatrixD invcov = Msys.Invert();

	chisqTest = 0;
	for(int i = 0; i < nBins; i++){
		for(int j = 0; j < nBins; j++){
			// (obsi-predi)*(invcov)*(obsj-predj)
			chisqTest += (testSpec[i] - predSpec.collapsed_vector[i])*invcov[i][j]*(testSpec[j] - predSpec.collapsed_vector[j]);
			// chisqTest += (testSpec.collapsed_vector[i] - predSpec.collapsed_vector[i])*invcov[i][j]*(testSpec.collapsed_vector[i] - predSpec.collapsed_vector[j]);

		}
	}
	return chisqTest;
}//end of GetChiSqFromSpectra

// int testfcn(std::vector<float> obsSpec,TMatrixD Msys,
// 				std::vector<std::tuple<double,double,double,int,int,int>> eventtuple_list,
// 				float m41, float ue4, float um4, bool debug = false){

void testfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

	// function to calculate the chi2 between fake universe and the mc events with osc weights
	float e_app = 4*pow(par[1],2)*pow(par[2],2);
	float e_dis = 4*pow(par[1],2)*(1-pow(par[1],2));
	float m_dis = 4*pow(par[2],2)*(1-pow(par[2],2));

	// find the closest mass spectra
	int lowidx, highidx;
	double prevval = 0;
	for(int i = 1;i<a_sinsqSpec.size();i++){
		// std::cout<<std::get<1>(a_sinsqSpec.at(i))<<" "<<par[0]<<" "<<prevval<<std::endl;
		if (par[0]==std::get<1>(a_sinsqSpec.at(i))){
			lowidx = i;
			highidx =i;
			break;
		}
		else if (par[0]<std::get<1>(a_sinsqSpec.at(i)) && par[0] >prevval ){
			lowidx = i-1;
			highidx =i;
			break;
		}
		else if( i == a_sinsqSpec.size()-1){
			lowidx = i;
			highidx =i;
		}
		else prevval = std::get<1>(a_sinsqSpec.at(i));
	}

	int closeidx = lowidx;
	if (lowidx <0) lowidx =0;
	if (lowidx<a_sinsqSpec.size()-2){
		double diff1 = par[0] - std::get<1>(a_sinsqSpec.at(lowidx));
		double diff2 = std::get<1>(a_sinsqSpec.at(highidx)) - par[0];
		if (diff2 <diff1) closeidx =highidx;
	}
	SBNspec inspec = std::get<0>(a_sinsqSpec.at(closeidx));
	SBNspec newSpec = GetOscillatedSpectra(cvSpec, inspec,e_app,e_dis, m_dis);

	cov = GetTotalCov(fakeData, newSpec, *covFracSys_collapsed);
	f =  GetChiSqFromSpectra(fakeData, newSpec, cov);
} // end of fcn function
