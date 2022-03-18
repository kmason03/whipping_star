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
SBNspec GetOscillatedSpectra(const SBNspec& cvSpec, SBNspec massSpec,
			float e_app, float e_dis, float m_dis);
TMatrixD GetTotalCov(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Mfracsys);
float GetLLHFromVector(const std::vector<float>& obsSpec, const SBNspec& expSpec,const TMatrixD& Msys,bool prints);
void testfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
std::string ZeroPadNumber(int num, int digits=3);

// define some global variables
std::string xml = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/xml/TotalThreePlusOne_full.xml";
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
SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_Bkg.SBNspec.root",xml);
std::array<double,nBins>  a_pred;
std::ofstream chifile;
std::vector<std::tuple<SBNspec,double>> a_sinsqSpec;
// TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/build/Appearence/DL_full.SBNcovar_bigbin.root","read");
TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins.SBNcovar.root","read");

TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");

int main(int argc, char* argv[]){
  // get the input integer: require 1
  int specific_entry = atoi(argv[1]);
	if (specific_entry > 15624) return 0;
	else{
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
		std::string textid = ZeroPadNumber(specific_entry, 4);
		chifile.open("chis_"+textid+".txt", std::ios_base::app);
		TFile *fout=new TFile("DLFCTests_hists.root","RECREATE");

		// Print binedges for easier plotting
		// if(printbins) printbinedges();

		// Load up the necesary bits and store them in vectors on the stack **
		cvSpec.Scale("fullosc",0.0);
		cvSpec.RemoveMCError();

		// make a tuple of the spectra and mass term
		for(int mi = 0; mi < 400; mi++){
			mnu = pow(10.,((mi+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			std::stringstream stream;
			stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
			std::cout<<mi<<std::endl;
			std::string infile = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_SINSQ_dm_"+stream.str()+".SBNspec.root";
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

		// get the oscillated spectra
		oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);


		 // I'll be generating a new universe around the oscillated spectrum
		// SBNchi TrueChi(cvSpec, true);
		SBNchi TrueChi(oscSpec, *covFracSys);
		fsys->cd();
		TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
		int b = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);


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
			TMatrixD cov_pT= GetTotalCov(fakeData, oscSpec, *oscFracSys_collapsed);
			float chi_pT = GetLLHFromVector(fakeData, oscSpec, cov_pT,false);

			// the minimizer from throw point
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
			Double_t dm2_start = mnu_base;
			Double_t ue4_start = ue_base;
			Double_t umu4_start = um_base;
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
			chifile<<chi_pT-fmin<<std::endl;
			delete gMinuit_simple;



		} //end of universeloop
		delete oscFracSys_collapsed;
	// chifile<< std::endl;
		return 0;
	}
} // end of main function


//-----------------------------------------------------------------------------
//----------------------HELPER FUNCTIONS---------------------------------------
//----------------------------------------------------------------------------

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
