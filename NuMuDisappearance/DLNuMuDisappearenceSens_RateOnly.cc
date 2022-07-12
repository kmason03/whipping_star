#include <iostream>
#include <sstream>
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

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "MillsFunctions.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"
#include "SBNcovariance.h"
#include "prob.h"
#include <assert.h>

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

// Globals:
// define some global variables
std::string xml = "/cluster/tufts/wongjiradlab/jmills09/whipping_star/xml/numudisappearance.xml";
std::string tag = "DL";
// set these parameters at the very start
const double dm2_lowbound(0.01), dm2_hibound(100);
const int dm2_grdpts(25);
double  mnu;


int main(int argc, char* argv[]){
	// testFunction();

	int iarg = 0;
	opterr=1;
	int index;
	bool gen = false;
	const double scaleFactor(1.0);
	// const double scaleFactor(2.0);
	// const double scaleFactor((190454.0/4848.0));

	bool numudis = false;
	bool nueapp = true;
	bool combined = false;
	int mass_start = -1;
	bool shapeonly = false;
	int miToDo = 0;
	bool doOneMass = false;
	bool statsOnly = false;
	const struct option longopts[] ={
	{"xml", 		required_argument, 	0, 'x'},
	{"gen",	no_argument, 0, 'g'},
	{"dis",	no_argument,0,'d'},
	{"app", no_argument,0,'a'},
	{"comb", no_argument,0,'c'},
	{"statsOnly", no_argument,0,'st'},
	{"part", required_argument,0,'p'},
	{"mibase", required_argument,0,'mi'},

	{0,			no_argument, 		0,  0},
	};

	while(iarg != -1){
		iarg = getopt_long(argc,argv, "x:dscp:g", longopts, &index);
		switch(iarg){
			case 'x':
				xml = optarg;
				break;
			case 'g':
				gen = true;
				break;
			case 'd':
				numudis = true;
				break;
			case 'a':
				nueapp = true;
				break;
			case 'c':
				combined = true;
				break;
			case 'p':
				mass_start = atoi(optarg);
				break;
			case 'mi':
				doOneMass = true;
				miToDo = atoi(optarg);
				break;
			case 'st':
				statsOnly = true;
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
			return 0;
		}
	}

	std::string tag = "DL";
	// set these parameters at the very start
	const double sin22th_lowbound(1e-2), sin22th_hibound(1.0);
	const int sin22th_grdpts(25);
	const int nFakeExp(1000);
	const int nBins(19);
	const int nBins_final(1);
	double  mnu, sin22th, ue, um,_mnu, _sin22th, _ue, _um, chisqTest, chisqMin, chisqPT;
	int count;
	std::vector<float> fakeData;
	fakeData.resize(nBins);
	SBNspec truespec_first, appSpec_first, truespec_second, appSpec_second, testSpec;
	std::array<double,nBins>  a_pred;

	//PART 1: precompute all of our sin and sin2 amplitudes so we don't need to later
	if(gen){
		// use true if we need to regenerate a covariance matrix
		if(false){
		    SBNcovariance _covar(xml);
		   _covar.FormCovarianceMatrix(tag);
		   _covar.PrintMatricies(tag);
		   _covar.frac_covariance.Print();
		}
		// for(int mi = 13; mi < 25; mi++){
		for(int mi = 0; mi < 13; mi++){
			generate_spectra(mi, xml, tag, dm2_lowbound, dm2_hibound, dm2_grdpts);
		}
		return 1;
	}

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	if(!gen){
		// save output to text files
		std::ofstream chifile;
		std::string chiFileName;
		if (statsOnly){
			if (!shapeonly){
				chiFileName = "rterms_rateonly_statsonly_1x_nonfreq.txt";
			}
			else{
				chiFileName = "rterms_rateonly_statsonly_shapeonly_1x_nonfreq.txt";
			}
			if (doOneMass){
				if (!shapeonly){
					chiFileName = "rterms_rateonly_statsonly_1x_nonfreq_pt_"+ZeroPadNumber(miToDo)+".txt";
				}
				else{
					chiFileName = "rterms_rateonly_statsonly_shapeonly_1x_nonfreq_pt_"+ZeroPadNumber(miToDo)+".txt";
				}
			}
		}
		else{
			if (!shapeonly){
				chiFileName = "rterms_rateonly_1x_nonfreq.txt";
			}
			else{
				chiFileName = "rterms_rateonly_shapeonly_1x_nonfreq.txt";
			}
			if (doOneMass){
				if (!shapeonly){
					chiFileName = "rterms_rateonly_1x_nonfreq_pt_"+ZeroPadNumber(miToDo)+".txt";
				}
				else{
					chiFileName = "rterms_rateonly_shapeonly_1x_nonfreq_pt_"+ZeroPadNumber(miToDo)+".txt";
				}
			}
		}

		chifile.open(chiFileName.c_str(), std::ios_base::app);
		dumpCoordFile(dm2_grdpts, dm2_lowbound, dm2_hibound, sin22th_grdpts, sin22th_lowbound);


		std::string endStr= "_Bkg.SBNspec.root";
		// if (scaleFactor == 1.0){
		// 	endStr = "_Bkg.SBNspec.root";
		// }
		// else{
		// 	endStr = "_1x_Bkg.SBNspec.root";
		// }
		SBNspec cvSpec(tag+endStr,xml);
		cvSpec.RemoveMCError();
		cvSpec.CollapseVector();
		//Scale all by ScaleFactor
		cvSpec.ScaleAll(scaleFactor);
		// cvSpec.Scale("fullosc",0.0);
		std::cout << "CV Spec:\n";
		cvSpec.PrintCollapsedVector();
		std::cout << "\n\n";
		// assert (1==2);

		std::vector<SBNspec> a_sinsqSpec;
		getPrecompSpec(a_sinsqSpec, dm2_grdpts, dm2_hibound, dm2_lowbound, xml, scaleFactor);
		//Bring in our covariance matrix!
		// Stats + sys
		// Load up cov matrix and add in detector variation component	**
		std::string covInFile;
		if (statsOnly) {covInFile = "JOSH_1m1p_StatOnly.SBNcovar.root";}
		else {covInFile = "JOSH_1m1p.SBNcovar.root";}
		TFile * fsys = new TFile(covInFile.c_str(),"read");
		TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
		TMatrixD * covFracSys_collapsed = (TMatrixD*)fsys->Get("frac_covariance_collapsed");

		std::cout<<"Loaded Covariance"<<std::endl;
		// // Load up oscillation model
		std::cout << "NUMU DISAPPEARANCE" << std::endl;
		float mnu_base, um_base, ue_base, sin22th_first;
		// first get the central grid point we are testing
		// Setting up a way to call the script piecemeal.
		int miBaseStart = 0;
		int miBaseLimit = dm2_grdpts;
		if (doOneMass){
			miBaseStart = miToDo;
			miBaseLimit = miToDo + 1;
		}
		double halfRTermMin = 9999999999;
		double chisqMin = 0;
		double detTermMin = 0;
		std::vector<double> detTermTest_v(dm2_grdpts*sin22th_grdpts,0);
		std::vector<double> chisqTest_v(dm2_grdpts*sin22th_grdpts,0);
		std::vector<double> halfRTermTest_v(dm2_grdpts*sin22th_grdpts,0);
		int vec_idx = -1;
		int min_mi = -1;
		int min_si = -1;
		for(int mi_base = miBaseStart; mi_base <miBaseLimit; mi_base++){
			for(int sin22thi_base = 0; sin22thi_base <sin22th_grdpts; sin22thi_base++){
				// if ((sin22thi_base != 24) || (mi_base != 15)) continue;
				vec_idx++;
				// std::cout << mi_base << " " << sin22thi_base << "\n";
				// if ((mi_base == 14) && (sin22thi_base == 0)){
				// 	assert (1==2);
				// }
				SBNspec thisGridCVSpec = cvSpec;
				thisGridCVSpec.RemoveMCError();
				TMatrixD cov(nBins,nBins), invcov(nBins_final,nBins_final), covCopy(nBins,nBins);
					//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
				sin22th_first = pow(10.,((sin22thi_base+0.5)/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
				mnu_base      = pow(10.,((mi_base+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
				// For NuMu Disappearance:
				um_base = 1.0 ; //THIS VAUE IS ALSO DUMMY
				ue_base = 1.0 ; // Dummy Value to satisfy code
				// current test model
				// std::cout << "NU Base Model: m41^2:" << pow(mnu_base,2) << " ue:" << ue_base << " sin2:" << sin22th_first << std::endl;

				// NuMu Disappearance Math described by equations 4.20-4.22 in Davio's thesis.
				truespec_first = cvSpec;
				// std::cout<<"cv spec"<<std::endl;
				// std::cout << "Before and After Spectra \n";
				// cvSpec.PrintCollapsedVector();

				scaleDisappearedSpec(truespec_first, a_sinsqSpec[mi_base], sin22th_first);
				// std::cout<<"oscillated spec"<<std::endl;
				// truespec_first.PrintCollapsedVector();

				SBNchi TrueChi(truespec_first,*covFracSys);
				TrueChi.ReloadCoreSpectrum(&truespec_first);
				TMatrixD covFracSys_base = getZeroedMat(nBins);
				TrueChi.FillCollapsedFractionalMatrix(&covFracSys_base);

				if(shapeonly){
					return -1; //Not configured for 1 bin test, that's the reason why we're doing this
					scaleCovShapeOnly(cov, covFracSys_base, thisGridCVSpec.collapsed_vector, truespec_first.collapsed_vector, nBins); // cov, covfraccollapsed, obs, exp(gets modified), nbins					// makeSpectraPlotsComp(cvSpec.collapsed_vector, truespec_first.collapsed_vector, sin22thi_base, mi_base, sin22th_first, pow(mnu_base,2));
				}
				// Shape+Rate
				else{
					scaleCovShapeRate(cov, covFracSys_base, thisGridCVSpec.collapsed_vector, truespec_first.collapsed_vector, nBins); // cov, covfraccollapsed, obs, exp(gets modified), nbins
				}
				// return 0;
				// covCopy = cov_single;
				// invcov = covCopy.Invert();

				std::vector<double> obs_v(1,0);
				std::vector<double> exp_v(1,0);

				for (int i=0;i<nBins;i++){
					obs_v[0] += thisGridCVSpec.collapsed_vector[i];
					exp_v[0] += truespec_first.collapsed_vector[i];
				}
				std::cout << obs_v[0] << " Observed\n";
				std::cout << exp_v[0] << " Expected\n";
				TMatrixD cov_single = getZeroedMat(nBins_final);
				for (int i=0;i<nBins;i++){
					for (int j=0;j<nBins;j++){
						cov_single[0][0] += cov[i][j];
					}
				}
				makeSpectraPlotsComp(obs_v, exp_v, sin22thi_base, mi_base, sin22th_first, pow(mnu_base,2), nBins_final, "onebin_expSpec");


				invcov[0][0] = 1.0/cov_single[0][0];
				std::cout << cov_single[0][0] << " Cov One Bin\n";
				std::cout << invcov[0][0] << " InvCov One Bin\n";

				double chisqTest;
				double detTermTest;
				double halfRTermTest;
				calcRTerm(detTermTest, chisqTest, obs_v, exp_v, invcov, cov_single, nBins_final); // obs, exp, inv, cov, nbins
				halfRTermTest = detTermTest+chisqTest;
				detTermTest_v[vec_idx] 		= detTermTest;
				chisqTest_v[vec_idx] 			= chisqTest;
				halfRTermTest_v[vec_idx] 	= halfRTermTest;

				if (((mi_base == 0) && (sin22thi_base ==0)) || ((mi_base == 23) && (sin22thi_base ==19))){
					std::cout << "mi_base, "<<mi_base << ", sin22, "<< sin22thi_base <<"\n";
					std::cout << detTermTest << ", " << chisqTest << ", " << detTermTest+chisqTest << "\n";
				}

				if (halfRTermTest < halfRTermMin){
					halfRTermMin = halfRTermTest;
					detTermMin = detTermTest;
					chisqMin = chisqTest;
					min_mi = mi_base;
					min_si = sin22thi_base;
				}
				if (halfRTermTest < 0){
					std::cout << halfRTermTest << " Is Negative\n";
					assert (1==2);
				}
			}// End of sin grid loop
		}// End of m grid loop
		double chisqTest;
		double detTermTest;
		double halfRTermTest;
		for (int ix=0;ix<detTermTest_v.size();ix++){
			chisqTest = chisqTest_v[ix];
			detTermTest = detTermTest_v[ix];
			halfRTermTest = halfRTermTest_v[ix];
			chifile << halfRTermTest-halfRTermMin << " ";
			std::cout << detTermTest << ", " << chisqTest << ", " << detTermMin << ", " << chisqMin << ", " << halfRTermTest-halfRTermMin << ", Term Breakdowns\n";
			chifile<< std::endl;
		}
		std::cout << "Minimum Match is:\n" << min_mi << " " << min_si <<"\n";

	}// End of not gen loop
	std::cout << "End of Script \n";
	return 0;
}
