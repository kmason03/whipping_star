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
SBNspec GetOscillatedSpectra(SBNspec cvSpec, SBNspec massSpec,
			float e_app, float e_dis, float m_dis);
TMatrixD GetTotalCov(SBNspec obsSpec, SBNspec expSpec, TMatrixD Mfracsys);
float GetLLHFromSpectra(SBNspec obsSpec, SBNspec expSpec, TMatrixD Msys);

// define some global variables
std::string xml = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/xml/TotalThreePlusOne_full.xml";
bool gen = false;
bool printbins = true;
int mass_start = -1;
std::string tag = "DL_full";
// set these parameters at the very start
const double dm2_lowbound(0.01), dm2_hibound(100);
const double ue4_lowbound(0.01), ue4_hibound(0.5);
const double umu4_lowbound(0.01), umu4_hibound(0.5);
// to genergate dm2_grdpts = 100
const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
const int nBins_e(12),nBins_mu(19);
const int nBins = nBins_e+nBins_mu;
TMatrixD cov_osc(nBins,nBins);
TMatrixD cov_cv(nBins,nBins);
double  mnu, ue, umu;
int count;
std::vector<float> fakeData;
SBNspec  appSpec, innerSpec, oscSpec;
SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/"+tag+".SBNspec.root",xml);
std::array<double,nBins>  a_pred;
std::ofstream coordfile;
std::ofstream chifile;;
Float_t z[5],x[5],y[5],errorz[5];
std::vector<std::tuple<SBNspec,double>> a_sinsqSpec;
TMatrixD * covFracSys_collapsed;

int main(){

	// open output text files
	coordfile.open("bins_sens_bigbin.txt", std::ios_base::app);
	chifile.open("chis_sens_bigbin.txt", std::ios_base::app);
	// Print binedges for easier plotting
	if(printbins) printbinedges();

	// Load up the necesary bits and store them in vectors on the stack **
	cvSpec.Scale("fullosc",0.0);
	cvSpec.CollapseVector();
	cvSpec.RemoveMCError();
	// cvSpec.PrintFullVector();
	// make a tuple of the spectra and mass term
	for(int mi = 0; mi < 400; mi++){
		mnu = pow(10.,((mi+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		std::stringstream stream;
		stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
		std::string infile = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_SINSQ_dm_"+stream.str()+".SBNspec.root";
		auto inspec = SBNspec(infile,xml);
		// inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
		inspec.CollapseVector();
		inspec.RemoveMCError();
		std::tuple<SBNspec,double> singletup (inspec,mnu);
		a_sinsqSpec.push_back(singletup);
	}

	// Stats + sys
	// Load up cov matrix and add in detector variation component	**
	// TFile * fsys = new TFile("katieversion_total.SBNcovar.root","read");
	// TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/build/Appearence/DL_full.SBNcovar_bigbin.root","read");
	TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");
	TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
	SBNchi cvChi(cvSpec, *covFracSys);
	TMatrixD * cvFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
	int a = cvChi.FillCollapsedFractionalMatrix(cvFracSys_collapsed);

	// first get the central grid point we are throwing universes around
	int idx=0;
	for(int mi_base = 18; mi_base <19; mi_base++){
		for(int uei_base = 0; uei_base < 1; uei_base++){
			for(int umui_base = 24; umui_base < 25; umui_base++){
				std::cout<<idx<<std::endl;
				idx+=1;
					//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
				int mi_base_new = mi_base*(400/dm2_grdpts);
				float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
				float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
				float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
				// calculate scaling factors
				float e_app = 4*pow(ue_base,2)*pow(um_base,2);
				float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
				float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));
				// current test model
				// std::cout << "NU Base Model: m41^2:" << mnu_base << " ue:" << ue_base << " um:" << um_base << std::endl;
				// std::cout <<"scale factors: "<<e_app<<" "<<e_dis<<" "<<m_dis<<std::endl;

				// get the oscillated spectra and cov matrix
				oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
				SBNchi TrueChi(cvSpec, *covFracSys);
				TrueChi.ReloadCoreSpectrum(&oscSpec);
				TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
				int run = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);

				// arguments: obs,exp,cov
				cov_osc= GetTotalCov(cvSpec, oscSpec, *oscFracSys_collapsed);
				// arguments: obs,exp,cov
				float llh_osc = GetLLHFromSpectra(cvSpec, oscSpec, cov_osc);
				//
				// cov_cv= GetTotalCov(cvSpec, cvSpec, *cvFracSys_collapsed);
				// float llh_cv = GetLLHFromSpectra(cvSpec, cvSpec, cov_cv);

				chifile<<llh_osc<<std::endl;
			}// end of loop over base umu
		}//end of loop over base ue
	}//end of loop over base mass

	return 0;
} // end of main function


//-----------------------------------------------------------------------------
//----------------------HELPER FUNCTIONS---------------------------------------
//----------------------------------------------------------------------------

void printbinedges(){
	// funtion that prints the bins to the output textfile
	for(int mi = 0; mi <= dm2_grdpts; mi++){
		mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		coordfile << pow(mnu,2) << " ";
	}
	coordfile << std::endl;

	for(int uei = 0; uei <= ue4_grdpts; uei++){
		ue = pow(10.,(uei/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
		coordfile << ue << " ";
	}
	coordfile << std::endl;

	for(int umui = 0; umui <= umu4_grdpts; umui++){
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
	massSpec.Scale("ext",0.0);
	// massSpec.PrintCollapsedVector();
	testSpec.Add(&massSpec);
	// std::cout<<"oscillated spec"<<std::endl;
	// testSpec.PrintCollapsedVector();
	return testSpec;
}//end of GetOscillatedSpectra

TMatrixD GetTotalCov(SBNspec obsSpec, SBNspec expSpec, TMatrixD Mfracsys){
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
					fullcov[i][j] += 3.0 / (1.0/obsSpec.collapsed_vector[i] + 2.0/expSpec.collapsed_vector[i]);
					// cov[i][j] += 3.0 / (1.0/testSpec.collapsed_vector[i] + 2.0/predSpec.collapsed_vector[i]);
				}
				else {
					fullcov[i][j] += expSpec.collapsed_vector[i]/2.0;
				}
			}
		}//end of first bin loop
	}//end of second bin loop
	return fullcov;
}//end of GetTotalCov

float GetLLHFromSpectra(SBNspec obsSpec, SBNspec expSpec, TMatrixD Msys){
	// // function to calculate a chi2 (shape + rate)
	// // inputs:
	// // obsSpec: "data" spectra
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
			chisqTest += (obsSpec.collapsed_vector[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec.collapsed_vector[j] - expSpec.collapsed_vector[j]);
			if(i==j)std::cout<<(obsSpec.collapsed_vector[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec.collapsed_vector[j] - expSpec.collapsed_vector[j])<<std::endl;
		}
	}
	std::cout<<chisqTest<<" "<<log(Msys.Determinant())<<std::endl;
	// now need ln(det(2Pi*M))
	// TMatrixD tempcov = 2*3.14159265358979323846*Msys;
	chisqTest += log(Msys.Determinant());

	return chisqTest;
}//end of GetLLHFromSpectr
