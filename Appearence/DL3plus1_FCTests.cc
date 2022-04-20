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
#include <gperftools/profiler.h>

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
void printbinedges();
SBNspec GetOscillatedSpectra(const SBNspec& cvSpec, SBNspec massSpec,
			float e_app, float e_dis, float m_dis);
TMatrixD GetTotalCov(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Mfracsys);
float GetLLHFromVector(const std::vector<float>& obsSpec, const SBNspec& expSpec,const TMatrixD& Msys,bool prints);
void testfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
std::string ZeroPadNumber(int num, int digits=3);

// define some global variables
std::string SBNHOME="/home/twongjirad/working/larbys/sbnfit/";
//std::string SBNHOME="/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/";
//std::string xml = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/xml/TotalThreePlusOne_full.xml";
std::string xml = SBNHOME+"/xml/TotalThreePlusOne_full_tmw.xml";  
bool draw = true;
bool printbins = true;
int mass_start = -1;
std::string tag = "DL_full";
// set some start parameters
const double dm2_lowbound(0.01), dm2_hibound(100);
const double ue4_lowbound(0.01), ue4_hibound(0.5);
const double umu4_lowbound(0.01), umu4_hibound(0.5);
// to genergate dm2_grdpts = 400
const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
const int nBins_e(12),nBins_mu(19);
const int nBins = nBins_e+nBins_mu;
const int nFakeExp(1);
double  mnu, ue, umu;
int count;
std::vector<float> fakeData;
TMatrixD cov(nBins,nBins);
SBNspec  appSpec, innerSpec, oscSpec;
// SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/"+tag+".SBNspec.root",xml);
//SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_Bkg.SBNspec.root",xml);
SBNspec cvSpec(SBNHOME+"/data/MassSpectraBigBins/"+tag+"_Bkg.SBNspec.root",xml);
std::array<double,nBins>  a_pred;
std::ofstream coordfile;
std::ofstream chifile;
std::ofstream covfile;
std::ofstream spacefile;
std::ofstream specfile;
std::vector<std::tuple<SBNspec,double>> a_sinsqSpec;
// TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/build/Appearence/DL_full.SBNcovar_bigbin.root","read");
// TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins.SBNcovar.root","read");
//TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");
TFile * fsys = new TFile((SBNHOME+"/data/systematics/katieversion_bigbins_tot.SBNcovar.root").c_str(),"read");

TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");

int main(int argc, char* argv[]){

    ProfilerStart("DL3plus1_FCTests.prof");
    
  // get the input integer: require 1
  int specific_entry = atoi(argv[1]);
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

	// resize fakeData vector
	fakeData.resize(nBins);

	// open output files
	std::string textid = ZeroPadNumber(specific_entry, 5);
	chifile.open("chis_"+textid+".txt", std::ios_base::app);
	covfile.open("cov_tot_big_220314.txt", std::ios_base::app);
	spacefile.open("space.txt", std::ios_base::app);
	specfile.open("spec_FC.txt", std::ios_base::app);
	TFile *fout=new TFile("DLFCTests_hists.root","RECREATE");

	// Print binedges for easier plotting
	// if(printbins) printbinedges();

	// Load up the necesary bits and store them in vectors on the stack **
	cvSpec.Scale("fullosc",0.0);
	// cvSpec.Scale("1e1p_nue",0.0);
	cvSpec.RemoveMCError();
	std::vector<double> cv_v = cvSpec.collapsed_vector;

	// make a tuple of the spectra and mass term
	for(int mi = 0; mi < 400; mi++){
		mnu = pow(10.,((mi+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		std::stringstream stream;
		stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
		//std::string infile = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_SINSQ_dm_"+stream.str()+".SBNspec.root";
		std::string infile = SBNHOME+"/data/MassSpectraBigBins/"+tag+"_SINSQ_dm_"+stream.str()+".SBNspec.root";
		std::cout<<"pre-calc dm2[" << mi << "]: " << infile << std::endl;		
		auto inspec = SBNspec(infile,xml);
		inspec.CollapseVector();
		inspec.RemoveMCError();
		std::tuple<SBNspec,double> singletup (inspec,mnu);
		a_sinsqSpec.push_back(singletup);
	}

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
	std::cout << "NU Base Model: m41:" <<pow(mnu_base,2) <<" m41^2:"<<mnu_base<< " ue:" << ue_base << " um:" << um_base << std::endl;
	std::cout<<"scale factors: "<<e_app<<" "<<e_dis<<" "<<m_dis<<std::endl;

	// get the oscillated spectra
	oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
	std::cout << "Full Vector of CV spectrum =========================" << std::endl;
	cvSpec.PrintFullVector(true);
	std::cout << "======================================================" << std::endl;	
	std::cout << "Matching dm2 parameter: " << std::get<1>(a_sinsqSpec.at(mi_base_new)) << std::endl;
	std::cout << "Full Vector of this spectrum =========================" << std::endl;
	std::get<0>(a_sinsqSpec.at(mi_base_new)).PrintFullVector(true);
	std::cout << "======================================================" << std::endl;
	std::vector<double> osc_v = oscSpec.collapsed_vector;
	// for(int i=0;i<nBins;i++){
	// 	specfile<<osc_v[i]<<" ";
	// }
	// specfile<<std::endl;

	if(draw){
		fout->cd();
		TH1D * oscSpec_1e1p_h = new TH1D("oscspec_1e1p","oscspec_1e1p ",nBins_e,0,nBins_e);
		TH1D * oscSpec_1mu1p_h = new TH1D("oscspec_1mu1p","oscspec_1mu1p",nBins_mu,0,nBins_mu);
		std::vector<double> oscspec_v=oscSpec.collapsed_vector;
		for(int i=0;i<nBins;i++){
			if (i<nBins_e) oscSpec_1e1p_h->SetBinContent(i,oscspec_v[i]);
			else oscSpec_1mu1p_h->SetBinContent(i-nBins_e,oscspec_v[i]);
		}

	}

	 // I'll be generating a new universe around the oscillated spectrum
	// SBNchi TrueChi(cvSpec, true);
	SBNchi TrueChi(oscSpec, *covFracSys);
	fsys->cd();
	TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("collapsed_frac_covariance");
	int b = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);
	for(short i = 0; i < nBins; i++){
		for(short j = 0; j < nBins; j++){
				covfile<< (*oscFracSys_collapsed)(i,j)<<" ";
		}
		covfile<<std::endl;
	}
	cvSpec.PrintCollapsedVector();

	oscSpec.PrintFullVector(true);

	// if ( true ) {
	//   return 0;
	// }
	
	TrueChi.pseudo_from_collapsed = true;
	TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomposition

	// Across several fake experiments for this grid point:
	for(int expi = 0; expi < nFakeExp; expi++){
		fakeData = TrueChi.GeneratePseudoExperiment();
		if(draw){
			fout->cd();
			TH1D * fakeSpec_1e1p_h = new TH1D("fakespec_1e1p","fakespec_1e1p ",nBins_e,0,nBins_e);
			TH1D * fakeSpec_1mu1p_h = new TH1D("fakespec_1mu1p","fakespec_1mu1p",nBins_mu,0,nBins_mu);
			for(int i=0;i<nBins;i++){
				if (i<nBins_e) fakeSpec_1e1p_h->SetBinContent(i,fakeData[i]);
				else fakeSpec_1mu1p_h->SetBinContent(i-nBins_e,fakeData[i]);
			}
		}


		// start of comptest
		// 1 get chi2pt
		float e_app_in = 4*pow(ue_base,2)*pow(um_base,2);
		float e_dis_in = 4*pow(ue_base,2)*(1-pow(ue_base,2));
		float m_dis_in = 4*pow(um_base,2)*(1-pow(um_base,2));
		TMatrixD cov_pT= GetTotalCov(fakeData, oscSpec, *oscFracSys_collapsed);
		float chi_pT = GetLLHFromVector(fakeData, oscSpec, cov_pT,true);
		chifile<<chi_pT<<" "<<mnu_base<<" "<<ue_base<<" "<<um_base<<std::endl;

		//2 get classic grid search
		float chi_min_grid =1000000;
		float m41_grid,ue_grid,umu4_grid;

		for(int mi_in = 0; mi_in <dm2_grdpts; mi_in++){
			std::cout<<mi_in<<std::endl;
			for(int uei_in = 0; uei_in < ue4_grdpts; uei_in++){
				for(int umui_in = 0; umui_in < umu4_grdpts; umui_in++){
					int mi_in_new = mi_in*(400/dm2_grdpts);
					float ue_val = pow(10.,(uei_in/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
					float um_val = pow(10.,(umui_in/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
					float mnu_val = pow(10.,((mi_in_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
					e_app_in = 4*pow(ue_val,2)*pow(um_val,2);
					e_dis_in = 4*pow(ue_val,2)*(1-pow(ue_val,2));
					m_dis_in = 4*pow(um_val,2)*(1-pow(um_val,2));
					// get the oscillated spectra
					SBNspec inSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_in_new)), e_app_in, e_dis_in, m_dis_in);
					SBNchi innerChi(inSpec, *covFracSys);
					TMatrixD * inFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
					int c = innerChi.FillCollapsedFractionalMatrix(inFracSys_collapsed);
					TMatrixD cov_grid = GetTotalCov(fakeData, inSpec, *inFracSys_collapsed);
					float chi_test = GetLLHFromVector(fakeData, inSpec, cov_grid, false);
					delete inFracSys_collapsed;
					// spacefile<<chi_test<<std::endl;
					// if(mi_in ==0 && uei_in==0 && umui_in==0){
					// 	for(int bin=0;bin<nBins;bin++){
					// 		specfile<<fakeData[bin]<<" ";
					// 	}
					// 	specfile<<std::endl;
					// }
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



		//3. just the minimizer from throw point
		// first some common variables
		Int_t ierflg = 0;
		Double_t dm2_step = .01;
		TString dm2_name = "m41";
		Double_t ue4_step = 0.001;
		TString ue4_name = "ue4";
		Double_t umu4_step = 0.001;
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
		arglist[0] = 500;
		arglist[1] = 3;

		// next initialize the minimizer
		TMinuit * gMinuit_simple= new TMinuit(3);
		gMinuit_simple->SetFCN(testfcn);
		Double_t dm2_start = m41_grid;
		Double_t ue4_start = ue_grid;
		Double_t umu4_start = umu4_grid;
		// start values near lower bound, but some wiggle room
		// initialize the parameters
		gMinuit_simple->mnparm(0,dm2_name,dm2_start,dm2_step,sqrt(dm2_lowbound),sqrt(dm2_hibound),ierflg);
		gMinuit_simple->mnparm(1,ue4_name,ue4_start,ue4_step,ue4_lowbound,ue4_hibound,ierflg);
		gMinuit_simple->mnparm(2,umu4_name,umu4_start,umu4_step,umu4_lowbound,umu4_hibound,ierflg);
		// run the minimizer
		gMinuit_simple->mnexcm("SEEk", arglist ,2,ierflg);
		// get the parameter outputs
		Double_t fmin,fedm,errdef;
		gMinuit_simple->mnpout(0,dm2_name,dm2_val,dm2_err,dm_low,dm_high,dm2_num);
		gMinuit_simple->mnpout(1,ue4_name,ue4_val,ue4_err,ue4_low,ue4_high,ue4_num);
		gMinuit_simple->mnpout(2,umu4_name,umu4_val,umu4_err,umu4_low,umu4_high,umu4_num);
		gMinuit_simple->mnstat(fmin,fedm,errdef,npari,nparx,istat);
		// std::cout<<"pt, fmin, coarse grid, fine grid"<<std::endl;
		// std::cout<<chi_pT<<" "<<fmin<<" "<<mnu_base<<" "<<ue_base<<" "<<um_base<<std::endl;
		chifile<<fmin<<" "<<dm2_val<<" "<<ue4_val<<" "<<umu4_val<<std::endl;

		if(draw){
			float e_app_m = 4*pow(ue4_val,2)*pow(umu4_val,2);
			float e_dis_m = 4*pow(ue4_val,2)*(1-pow(ue4_val,2));
			float m_dis_m = 4*pow(umu4_val,2)*(1-pow(umu4_val,2));

			// find the closest mass spectra
			int lowidx, highidx;
			double prevval = 0;
			for(int i = 1;i<a_sinsqSpec.size();i++){
				// std::cout<<std::get<1>(a_sinsqSpec.at(i))<<" "<<par[0]<<" "<<prevval<<std::endl;
				if (dm2_val==std::get<1>(a_sinsqSpec.at(i))){
					lowidx = i;
					highidx =i;
					break;
				}
				else if (dm2_val <std::get<1>(a_sinsqSpec.at(i)) && dm2_val >prevval ){
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
				double diff1 = dm2_val - std::get<1>(a_sinsqSpec.at(lowidx));
				double diff2 = std::get<1>(a_sinsqSpec.at(highidx)) - dm2_val;
				if (diff2 <diff1) closeidx =highidx;
			}
			// get spectra and covar
			SBNspec inspec = std::get<0>(a_sinsqSpec.at(closeidx));
			SBNspec newSpec = GetOscillatedSpectra(cvSpec, inspec,e_app_m,e_dis_m, m_dis_m);
			fout->cd();
			TH1D * bestSpec_1e1p_h = new TH1D("bestspec_1e1p","bestspec_1e1p ",nBins_e,0,nBins_e);
			TH1D * bestSpec_1mu1p_h = new TH1D("bestspec_1mu1p","bestspec_1mu1p",nBins_mu,0,nBins_mu);
			std::vector<double> bestspec_v=newSpec.collapsed_vector;
			for(int i=0;i<nBins;i++){
				if (i<nBins_e) bestSpec_1e1p_h->SetBinContent(i,bestspec_v[i]);
				else bestSpec_1mu1p_h->SetBinContent(i-nBins_e,bestspec_v[i]);
			}
			SBNchi tmpChi(newSpec, *covFracSys);
			TMatrixD * tmpFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
			int b = tmpChi.FillCollapsedFractionalMatrix(tmpFracSys_collapsed);
			TMatrixD tmpcov = GetTotalCov(fakeData, newSpec, *tmpFracSys_collapsed);
			// calculate -2LLH
			float testchi=  GetLLHFromVector(fakeData, newSpec, tmpcov, true);
			delete tmpFracSys_collapsed;
		}

		delete gMinuit_simple;



	} //end of universeloop
	std::cout<<mnu_base<<" "<<ue_base<<" "<<um_base<<std::endl;
	fout->Write();
	fout->Close();
	delete oscFracSys_collapsed;

	ProfilerStop();
	
// chifile<< std::endl;
	return 0;
} // end of main function


//-----------------------------------------------------------------------------
//----------------------HELPER FUNCTIONS---------------------------------------
//----------------------------------------------------------------------------

void printbinedges(){
	// funtion that prints the bins to the output textfile
	for(int mi = 0; mi < 400; mi++){
		mnu = pow(10.,((mi+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		coordfile << pow(mnu,2) << " ";
	}
	coordfile << std::endl;

	for(int uei = 0; uei <400; uei++){
		ue = pow(10.,(uei/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
		coordfile << ue << " ";
	}
	coordfile << std::endl;

	for(int umui = 0; umui < 400; umui++){
		umu = pow(10.,(umui/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
		coordfile << umu << " ";
	}
	coordfile << std::endl;
	return;
}//end of print bins function


SBNspec GetOscillatedSpectra(const SBNspec& cvSpec, SBNspec massSpec,
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
	massSpec.Scale("ext",0.0);
	testSpec.Add(&massSpec);

	return testSpec;
}//end of GetOscillatedSpectra

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
			if(i==j){
				if (expSpec.collapsed_vector[i] >0 ){
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
	float chisqTest;

	// inv cov for chi2calc
	TMatrixD invcov = Msys;
	invcov.Invert();

	// add the chi2-like part
	chisqTest = 0;
	for(int i = 0; i < nBins; i++){
		for(int j = 0; j < nBins; j++){
			// (obsi-predi)*(invcov)*(obsj-predj)
			if(i==j && prints) std::cout<<i<<" "<<obsSpec[i]<<" "<<expSpec.collapsed_vector[i]<<" "<<Msys[i][j]<<" "<<((obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]))<<std::endl;
			chisqTest += (obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]);
		}
	}
	// now need ln(det(2Pi*M))
	// TMatrixD tempcov = 2*3.14159265358979323846*Msys;
	// std::cout<<"chi2: "<<chisqTest<<" det: "<<log(tempcov.Determinant())<<std::endl;
	chisqTest += log(Msys.Determinant());

	return chisqTest;
}//end of GetLLHFromSpectra

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
	// get spectra and covar
	SBNspec inspec = std::get<0>(a_sinsqSpec.at(closeidx));
	SBNspec newSpec = GetOscillatedSpectra(cvSpec, inspec,e_app,e_dis, m_dis);
	SBNchi tmpChi(newSpec, *covFracSys);
	TMatrixD * tmpFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
	int b = tmpChi.FillCollapsedFractionalMatrix(tmpFracSys_collapsed);
	TMatrixD tmpcov = GetTotalCov(fakeData, newSpec, *tmpFracSys_collapsed);
	// calculate -2LLH
	f =  GetLLHFromVector(fakeData, newSpec, tmpcov, false);
	delete tmpFracSys_collapsed;
} // end of fcn function

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
