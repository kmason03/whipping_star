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
	std::string xml = "/uboone/app/users/kmason/whipping_star/xml/nueappearence.xml";
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
	const double sin22th_lowbound(1e-4), sin22th_hibound(1.0);
	const int dm2_grdpts(25), sin22th_grdpts(25);
	const int nFakeExp(1000);
	const int nBins(22);
	double  mnu, sin22th, ue, um,_mnu, _sin22th, _ue, _um, chisqTest, chisqMin, chisqPT;
	int count;
	TMatrixD cov(nBins,nBins), invcov(nBins,nBins);
	std::vector<float> fakeData;
	fakeData.resize(nBins);
	SBNspec trueSpec, appSpec, innerSpec;
	std::array<double,nBins>  a_pred;

	std::string s_fracDetSysMatrix = "/uboone/app/users/kmason/whipping_star/build/Appearence/fracdetvar.txt";
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
			mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

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
		return 1;
	}

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	if(!gen){
		// save output to text files
		std::ofstream coordfile;
		std::ofstream chifile;
		std::ofstream specfile_low;
		std::ofstream specfile_high;

		coordfile.open("bins.txt", std::ios_base::app);//std::ios_base::app
		chifile.open("chis.txt", std::ios_base::app);
		specfile_low.open("specs_low.txt", std::ios_base::app);
		specfile_high.open("specs_high.txt", std::ios_base::app);
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
		cvSpec.Scale("fullosc",0.0);
		cvSpec.CollapseVector();
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
		// mcstat error - found in selection grid from python (1/sqrt(n)) for each bin)
		std::vector<float> stkerr={0.082, 0.039, 0.031, 0.028, 0.028, 0.029, 0.031, 0.037, 0.045, 0.058, 0.061, 0.071, 0.083, 0.093, 0.137, 0.145, 0.160, 0.204, 0.316, 0.377, 0.408, 1.0};


		// // Load up oscillation model
		// lets just look at nue app for now

		// If we're doing nue appearance
		std::cout << "NUE APPEARANCE" << std::endl;
		float mnu_base, um_base, ue_base, sin22th_base;
		// first get the central grid point we are throwing universes around
		for(int mi_base = 0; mi_base <sin22th_grdpts; mi_base++){
			for(int sin22thi_base = 0; sin22thi_base <sin22th_grdpts; sin22thi_base++){
		// for(int mi_base = 0; mi_base <1; mi_base++){
		// 	for(int sin22thi_base = 23; sin22thi_base <24; sin22thi_base++){
					//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
				sin22th_base = pow(10.,(sin22thi_base/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
				mnu_base = pow(10.,(mi_base/float(dm2_grdpts)*TMath::Log10(10./.1) + TMath::Log10(.1)));
				um_base = .135;	// global best fit
				ue_base = pow(sin22th_base/float(4*um_base*um_base),.5);
				// current test model
				std::cout << "NU Base Model: m41^2:" << pow(mnu_base,2) << " ue:" << ue_base << " sin2:" << sin22th_base << std::endl;

				trueSpec = cvSpec;
				trueSpec.Scale("fullosc",0.0);
				std::cout<<"cv spec"<<std::endl;
				trueSpec.PrintCollapsedVector();
				appSpec = a_sinsqSpec[mi_base];
				appSpec.Scale("fullosc",sin22th_base);
				appSpec.Scale("bnb",0.0);
				appSpec.Scale("nue",0.0);
				appSpec.Scale("ccpi0",0.0);
				appSpec.Scale("ncpi0",0.0);
				appSpec.Scale("ext",0.0);

				trueSpec.Add(&appSpec);
				std::cout<<"oscillated spec"<<std::endl;
				trueSpec.PrintCollapsedVector();

				// for(int i = 0; i < nBins; i++){
				// 	specfile_low<<trueSpec.collapsed_vector[i]<<" ";
				// 	specfile_high<<trueSpec.collapsed_vector[i]<<" ";
				// }
				// specfile_high<<"\n";
				// specfile_low<<"\n";

				//  I'll be generating a new universe around the oscillated spectrum
				SBNchi TrueChi(trueSpec,*covFracSys);
				TrueChi.ReloadCoreSpectrum(&trueSpec);
				TrueChi.pseudo_from_collapsed = true;
				TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomposition


				// Across several fake experiments for this grid point:
				for(int expi = 0; expi < nFakeExp; expi++){
					fakeData = TrueChi.GeneratePseudoExperiment();	// [fakedata] = vector of floats

					if (expi==0){
						std::cout<<"osc spec"<<std::endl;
						trueSpec.PrintCollapsedVector();
					}

					// now compare this univers to every grid point to find the minimum chi2. Thats what we save
					// initialize as a really large number
					float minchi = 99999999;
					float chi_pT = 99999999;
					for(int mi_inner = 0; mi_inner <sin22th_grdpts; mi_inner++){

						// reset cov to zeros so it resets
						for(int i = 0; i < nBins; i++){
							for(int j = 0; j < nBins; j++){
								cov[i][j]=0.0;
							}
						}

						for(int sin22th_inner = 0; sin22th_inner <sin22th_grdpts; sin22th_inner++){
							// get the right scaling
							double sin22th_innerbase = pow(10.,(sin22th_inner/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));

							// get the prediction for the new grid point
							innerSpec =  a_sinsqSpec[mi_inner];
							innerSpec.Scale("fullosc",sin22th_innerbase);
							innerSpec.Scale("bnb",0.0);
							innerSpec.Scale("nue",0.0);
							innerSpec.Scale("ccpi0",0.0);
							innerSpec.Scale("ncpi0",0.0);
							innerSpec.Scale("ext",0.0);
							innerSpec.Add(&cvSpec);
							// innerSpec.PrintCollapsedVector();


							// calculate the full covariance matrix
							for(int i = 0; i < nBins; i++){
								for(int j = 0; j < nBins; j++){
										// remove some nans
										if(innerSpec.collapsed_vector[i] > 0 && innerSpec.collapsed_vector[j] >0 )
										cov[i][j] += (*covFracSys_collapsed)[i][j]*innerSpec.collapsed_vector[i]*innerSpec.collapsed_vector[j];

										if(i==j){
											// cov[i][j] += trueSpec.collapsed_vector[i];
											// add in stat errors start with CNP for "data" errors
											if (innerSpec.collapsed_vector[i] >0 ){
												cov[i][j] += 3.0 / (1.0/fakeData[i] + 2.0/innerSpec.collapsed_vector[i]);
												// poisson version
												// cov = (M-mu)**2 / (2*(mu - M + M*np.log(M/mu)))
												// cov[i][j] += pow(fakeData[i]-trueSpec.collapsed_vector[i],2) / (2.0*(trueSpec.collapsed_vector[i]- fakeData[i]+ fakeData[i]*log(fakeData[i]/trueSpec.collapsed_vector[i])));

											}
											else {
												cov[i][j] += innerSpec.collapsed_vector[i]/2.0;
											}
											// mc stat error
											cov[i][j] +=stkerr[i]*stkerr[i];
										}
									}
							}

							// inv cov for chi2calc
							invcov = cov.Invert();

							chisqTest = 0;
							for(int i = 0; i < nBins; i++){
								for(int j = 0; j < nBins; j++){
									// (obsi-predi)*(invcov)*(obsj-predj)
									chisqTest += (fakeData[i] - innerSpec.collapsed_vector[i])*invcov[i][j]*(fakeData[j] - innerSpec.collapsed_vector[j]);
								}
							}
							// std::cout<<chisqTest<<std::endl;
							// if its the pt the universe came from, save as pT
							if (mi_base == mi_inner && sin22thi_base ==sin22th_inner  ){
								chi_pT=chisqTest;
							}
							// is it the new minimum?
							if (chisqTest < minchi){
								minchi=chisqTest;
							}
						}//end of inner sin Loop
					}//end of inner mass loop


				// writing some test spectra to file for debugging
				// if (chisqTest<40){
				// 	for(int i = 0; i < nBins; i++){
				// 		specfile_low<<fakeData[i]<<" ";
				// 	}
				// 	specfile_low<<"\n";
				// }
				// if (chisqTest>70){
				// 	for(int i = 0; i < nBins; i++){
				// 		specfile_high<<fakeData[i]<<" ";
				// 	}
				// 	specfile_high<<"\n";
				// }

				// std::cout<<"chi2: "<<chisqTest<<std::endl;
				chifile << (chi_pT-minchi)<< " ";
			} //end of universeloop
			chifile<< std::endl;


			}
		}
	}
	return 0;
}
