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

int main(int argc, char* argv[]){
	std::string xml = "/uboone/app/users/jmills/whipping_star/xml/numudisappearance.xml";
	int iarg = 0;
	opterr=1;
	int index;
	bool gen = false;
	bool numudis = false;
	bool nueapp = true;
	bool combined = false;
	int mass_start = -1;
	bool shapeonly = false;
	const struct option longopts[] ={
	{"xml", 		required_argument, 	0, 'x'},
	{"gen",	no_argument, 0, 'g'},
	{"dis",	no_argument,0,'d'},
	{"app", no_argument,0,'a'},
	{"comb", no_argument,0,'c'},
	{"part", required_argument,0,'p'},
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
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
			return 0;
		}
	}

	std::string tag = "DL";
	// set these parameters at the very start
	const double dm2_lowbound(0.01), dm2_hibound(100);
	const double sin22th_lowbound(1e-2), sin22th_hibound(1.0);
	const int dm2_grdpts(25), sin22th_grdpts(25);
	const int nFakeExp(1);
	const int nBins(19);
	double  mnu, sin22th, ue, um,_mnu, _sin22th, _ue, _um, chisqTest, chisqMin, chisqPT;
	int count;
	TMatrixD cov(nBins,nBins), invcov(nBins,nBins);
	std::vector<float> fakeData;
	fakeData.resize(nBins);
	SBNspec trueSpec, appSpec, testSpec;
	std::array<double,nBins>  a_pred;

	std::string s_fracDetSysMatrix = "/uboone/app/users/jmills/whipping_star/NuMuDisappearance/fracdetvar_19bins.txt";
	//PART 1: precompute all of our sin and sin2 amplitudes so we don't need to later
	if(gen){
		// use true if we need to regenerate a covariance matrix
		if(true){
		    SBNcovariance _covar(xml);
		   _covar.FormCovarianceMatrix(tag);
		   _covar.PrintMatricies(tag);
		   _covar.frac_covariance.Print();
		}
		// start with a null model
		NeutrinoModel nullModel(0, 0, 0);
		SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
		SBNspec bkg = bkgo->spec_central_value;
		//SBNspec bkg("DL.SBNspec.root",xml);
		// set full osc to 0 for bkg spectrum
		// bkg.Scale("fullosc",0.0);
		bkg.WriteOut(tag+"_Bkg");

		float mnu;
		for(int mi = 13; mi < 25; mi++){
		// for(int mi = 0; mi < 13; mi++){

			mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

			//Model: mnu, ue4, um4
			//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
			//see davio's thesis eqs 4.14-16
			//we're precomputing, so we don't really care about the u's set to sin2=1
			// m41, Ue, Um
			// NeutrinoModel testModel(mnu, sqrt(.5), sqrt(.5)); //Nue App Maximized
			NeutrinoModel testModel(mnu, sqrt(.5), sqrt(.5)); //NuMu Disap Maximized

			// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
			SBNgenerate * gen = new SBNgenerate(xml,testModel);

			// Write them to files
			gen->WritePrecomputedOscSpecs(tag);
		}
		return 1;
	}

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	if(!gen){
		// save output to text files
		std::ofstream coordfile;
		std::ofstream chifile;

		coordfile.open("bins.txt", std::ios_base::app);//std::ios_base::app
		chifile.open("chis.txt", std::ios_base::app);
		// Print binedges for easier plotting **
		for(int mi = 0; mi <= dm2_grdpts; mi++){
			mnu = pow(10.,((mi)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			coordfile << pow(mnu,2) << " ";
		}
		coordfile << std::endl;

		for(int si = 0; si <= sin22th_grdpts; si++){
			sin22th = pow(10.,(si/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
			coordfile << sin22th << " ";
		}
		coordfile << std::endl;

		// Load up the necesary bits and store them in vectors on the stack **
		SBNspec cvSpec(tag+"_Bkg.SBNspec.root",xml);
		cvSpec.CollapseVector();
		// cvSpec.Scale("fullosc",0.0);
		// cvSpec.PrintFullVector();
		std::array<SBNspec,dm2_grdpts> a_sinsqSpec;
		for(int mi = 0; mi < dm2_grdpts; mi++){
			mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			std::stringstream stream;
			stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
			std::string infile = "DL_SINSQ_dm_"+stream.str()+".SBNspec.root";
			auto inspec = SBNspec(infile,xml);
			// inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
			inspec.CollapseVector();
			a_sinsqSpec[mi] = inspec;
		}

		//Bring in our covariance matrix!
		// Stats only
		SBNchi uboone_chi_statsonly(cvSpec,true);

		// Stats + sys
		// Load up cov matrix and add in detector variation component	**
		TFile * fsys = new TFile("DL.SBNcovar.root","read");
		TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
		TMatrixD * covFracSys_collapsed = (TMatrixD*)fsys->Get("collapsed_frac_covariance");
		// need to add in detvar still...
		std::ifstream file;
		double ddummy;
		file.open(s_fracDetSysMatrix);
		if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << s_fracDetSysMatrix << std::endl;
		for(short i = 0; i < nBins; i++){
			for(short j = 0; j < nBins; j++){
				file >> ddummy;
				(*covFracSys_collapsed)(i,j) += ddummy;
			}
		}
		file.close();
		std::cout<<"Loaded Covariance"<<std::endl;

		// // Load up oscillation model
		// lets just look at nue app for now

		// If we're doing nue appearance
		std::cout << "NUE APPEARANCE" << std::endl;
		std::cout << "NUMU DISAPPEARANCE" << std::endl;

		float mnu_base, um_base, ue_base, sin22th_base;
		// first get the central grid point we are testing
		for(int mi_base = 0; mi_base <sin22th_grdpts; mi_base++){
			for(int sin22thi_base = 0; sin22thi_base <sin22th_grdpts; sin22thi_base++){
					//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
				sin22th_base = pow(10.,(sin22thi_base/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
				mnu_base = pow(10.,(mi_base/float(dm2_grdpts)*TMath::Log10(10./.1) + TMath::Log10(.1)));
				// For Nue Appearance:
				// um_base = .135;	// global best fit
				// ue_base = pow(sin22th_base/float(4*um_base*um_base),.5);
				// For NuMu Disappearance:
				um_base = 1.0 ; //THIS VAUE IS ALSO DUMMY
				ue_base = 1.0 ; // Dummy Value to satisfy code
				// current test model
				std::cout << "NU Base Model: m41^2:" << pow(mnu_base,2) << " ue:" << ue_base << " sin2:" << sin22th_base << std::endl;

				// NuMu Disappearance Math described by equations 4.20-4.22 in Davio's thesis.
				trueSpec = cvSpec;
				// Davio
				// trueSpec = cvSpec;
				// disSpec = a_sinsqSpec[pTdm2i];
				// disSpec.Scale("bnb",-sin22th);
				// trueSpec.Add(&disSpec);
				// trueSpec.CollapseVector();

				trueSpec.Scale("fullosc",0.0);
				std::cout<<"cv spec"<<std::endl;
				std::cout << "Before and After Spectra \n";
				trueSpec.PrintCollapsedVector();
				appSpec = a_sinsqSpec[mi_base];
				// appSpec.Scale("fullosc",sin22th_base);
				appSpec.Scale("bnb",-1.0*sin22th_base);
				appSpec.Scale("nue",0.0);
				appSpec.Scale("ccpi0",-1.0*sin22th_base);
				appSpec.Scale("ncpi0",-1.0*sin22th_base);
				appSpec.Scale("ext",0.0);
				trueSpec.Add(&appSpec);
				std::cout<<"oscillated spec"<<std::endl;
				trueSpec.PrintCollapsedVector();

				//  I'll be generating a new universe around the oscillated spectrum
				SBNchi TrueChi(trueSpec,*covFracSys);
				TrueChi.ReloadCoreSpectrum(&trueSpec);
				TrueChi.pseudo_from_collapsed = true;
				TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomposition



				// Across several fake experiments:
				for(int expi = 0; expi < nFakeExp; expi++){
					fakeData = TrueChi.GeneratePseudoExperiment();	// [fakedata] = vector of floats
					// remove zero bins, they cause nans in the chi2 calc
					for(int i = 0; i < nBins; i++){
						if (fakeData[i] ==0) fakeData[i]=0.00000001;
					}
					// print out new oscillated Spectrum
					std::cout<<"osc spec"<<std::endl;
					trueSpec.PrintCollapsedVector();
					std::cout<<"new uni spec: ";
					for(int i = 0; i < nBins; i++){
						std::cout<<fakeData[i]<<" ";
					}
					std::cout<<"\n";

					// mcstat error - found in selection grid from python, in PlottingScripts.py 1/np.sqrt(yerrsq_mc) for each bin)
					// First bin is inf in 20 bin enu analysis, set error to 0.4 manually
					// std::vector<float> stkerr={0.4, 0.39555943, 0.14108268, 0.12801097, 0.10042527, 0.11516246, 0.1107993, 0.10101455, 0.09678998, 0.11705153, 0.11300651, 0.11581347, 0.14006509, 0.13711131, 0.1510809, 0.17920406, 0.17735261, 0.24390278, 0.26385786, 0.27852506};
					std::vector<float> stkerr={0.39555943, 0.14108268, 0.12801097, 0.10042527, 0.11516246, 0.1107993, 0.10101455, 0.09678998, 0.11705153, 0.11300651, 0.11581347, 0.14006509, 0.13711131, 0.1510809, 0.17920406, 0.17735261, 0.24390278, 0.26385786, 0.27852506};

					// haven't touched shape only yet, but leaving in from davio for use later.
					if(shapeonly){
						double covIntegral(0.0), predIntegral(0.0), obsIntegral(0.0),fnorm;
						for(int i = 0; i < nBins; i++){
							obsIntegral+=a_pred[i];
							predIntegral+=fakeData[i];
							for(int j = 0; j < nBins; j++){
								covIntegral += (*covFracSys)[i][j]*cvSpec.collapsed_vector[i]*cvSpec.collapsed_vector[j];
							}
						}
						fnorm = covIntegral/pow(predIntegral,2);
						for(int i = 0; i < nBins; i++){
							a_pred[i] *= (obsIntegral/predIntegral); // normalize prediction
						}
						for(int i = 0; i < nBins; i++){
							for(int j = 0; j < nBins; j++){
								cov[i][j] = ((*covFracSys)[i][j]-fnorm)*a_pred[i]*a_pred[j];
								if(i==j)
									cov[i][j] += a_pred[i];
							}
						}
					}
					// Shape+Rate
					else{
						for(int i = 0; i < nBins; i++){
							for(int j = 0; j < nBins; j++){
									// remove some nans
									if(trueSpec.collapsed_vector[i] > 0 && trueSpec.collapsed_vector[j] >0 )
									cov[i][j] += (*covFracSys_collapsed)[i][j]*trueSpec.collapsed_vector[i]*trueSpec.collapsed_vector[j];

									if(i==j){
										// cov[i][j] += trueSpec.collapsed_vector[i]; // Davio includes this
										// add in stat errors start with CNP for "data" errors
										if (trueSpec.collapsed_vector[i] >0 && fakeData[i] >0){
											// cov[i][j] += 3.0 / (1.0/fakeData[i] + 2.0/trueSpec.collapsed_vector[i]);
											// poisson version
											// cov = (M-mu)**2 / (2*(mu - M + M*np.log(M/mu)))
											cov[i][j] += pow(fakeData[i]-trueSpec.collapsed_vector[i],2) / (2.0*(trueSpec.collapsed_vector[i]- fakeData[i]+ fakeData[i]*log(fakeData[i]/trueSpec.collapsed_vector[i])));

										}
										else {
											cov[i][j] += trueSpec.collapsed_vector[i]/2.0;
										}
										// mc stat error
										cov[i][j] +=stkerr[i]*stkerr[i];
									}
								}
						}
					} //end of shape only
					// inv cov for chi2calc
					invcov = cov.Invert();

					chisqTest = 0;
					for(int i = 0; i < nBins; i++){
						for(int j = 0; j < nBins; j++){
							// (obsi-predi)*(invcov)*(obsj-predj)
							chisqTest += (fakeData[i] - trueSpec.collapsed_vector[i])*invcov[i][j]*(fakeData[j] - trueSpec.collapsed_vector[j]);
						}
					}

				chifile << chisqTest<< " ";
			}
			chifile<< std::endl;


			}
		}
	}
	std::cout << "End of Script \n";
	return 0;
}
