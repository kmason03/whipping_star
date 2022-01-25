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
	std::string xml = "/uboone/app/users/kmason/whipping_star/xml/TotalThreePlusOne.xml";
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

	std::string tag = "DL_full";
	// set these parameters at the very start
	const double dm2_lowbound(0.01), dm2_hibound(100);
	const double ue4_lowbound(0.0707), ue4_hibound(0.5);
	const double umu4_lowbound(0.0707), umu4_hibound(0.5);
	const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
	const int nBins_e(10),nBins_mu(19);
	const int nBins = nBins_e+nBins_mu;
	double  mnu, ue, umu, chisqTest;
	int count;
	TMatrixD cov(nBins,nBins), invcov(nBins,nBins);
	std::vector<float> fakeData;
	fakeData.resize(nBins);
	SBNspec trueSpec, appSpec, innerSpec;
	std::array<double,nBins>  a_pred;

	std::string s_fracDetSysMatrix = "/uboone/app/users/kmason/whipping_star/build/Appearence/fracdetvar_LEE.txt";
	std::string s_fracDetSysMatrix_mu = "/uboone/app/users/jmills/whipping_star/NuMuDisappearance/fracdetvar_19bins.txt";
	//PART 1: precompute all of our sin and sin2 amplitudes so we don't need to later
	if(gen){
		// use true if we need to regenerate a covariance matrix
		if(false){
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
		std::ofstream covfile;
		std::ofstream covtotalfile;
		std::ofstream covinvfile;
		std::ofstream cvspecfile;
		std::ofstream oscspecfile;

		coordfile.open("bins_full.txt", std::ios_base::app);
		chifile.open("chis_full.txt", std::ios_base::app);
		covfile.open("cov_full.txt", std::ios_base::app);
		covtotalfile.open("covtotal_full.txt", std::ios_base::app);
		covinvfile.open("covinv_full.txt", std::ios_base::app);
		cvspecfile.open("cvspec_full.txt", std::ios_base::app);
		oscspecfile.open("oscspec_full.txt", std::ios_base::app);

		// Print binedges for easier plotting **
		for(int mi = 0; mi <= dm2_grdpts; mi++){
			mnu = pow(10.,((mi)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
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

		// Load up the necesary bits and store them in vectors on the stack **
		SBNspec cvSpec("MassSpectraData/Full/"+tag+"_Bkg.SBNspec.root",xml);
		cvSpec.Scale("1e1p_fullosc",0.0);

		std::array<SBNspec,dm2_grdpts> a_sinsqSpec;
		for(int mi = 0; mi < dm2_grdpts; mi++){
			mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			std::stringstream stream;
			stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
			std::string infile = "MassSpectraData/Full/"+tag+"_SINSQ_dm_"+stream.str()+".SBNspec.root";
			auto inspec = SBNspec(infile,xml);
			// inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
			inspec.CollapseVector();
			a_sinsqSpec[mi] = inspec;
		}

		//Bring in our covariance matrix!
		// Stats + sys
		// Load up cov matrix and add in detector variation component	**
		TFile * fsys = new TFile("katieversion_total.SBNcovar.root","read");
		// TFile * fsys = new TFile("h1_v48_total.SBNcovar.root","read");
		TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
		// switch to collapsed version
		TMatrixD * covFracSys_collapsed = (TMatrixD*)fsys->Get("frac_covariance_collapsed");

		//save covar for plotting purposes
		for(short i = 0; i < nBins; i++){
			for(short j = 0; j < nBins; j++){
					covfile<< (*covFracSys_collapsed)(i,j)<<" ";
			}
			covfile<<std::endl;
		}
		std::cout<<"Loaded Covariance"<<std::endl;

		cvSpec.CollapseVector();
		cvSpec.PrintCollapsedVector();
		for(int i = 0; i < nBins; i++){
			cvspecfile<<cvSpec.collapsed_vector[i]<<" ";
		}
		cvspecfile<<"\n";


		// // Load up oscillation model
		// look at all modes
		float mnu_base, um_base, ue_base, sin22th_base;
		// first get the central grid point we are throwing universes around
		for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
			for(int uei_base = 0; uei_base < ue4_grdpts; uei_base++){
				for(int umui_base = 0; umui_base < umu4_grdpts; umui_base++){
						//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
					ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
					um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
					mnu_base = pow(10.,(mi_base/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(.1)));
					// calculate scaling factors
					float e_app = 4*pow(ue_base,2)*pow(um_base,2);
					float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
					float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));
					// current test model
					std::cout << "NU Base Model: m41^2:" << mnu_base << " ue:" << ue_base << " um:" << um_base << std::endl;
					std::cout<<"scale factors: "<<e_app<<" "<<e_dis<<" "<<m_dis<<std::endl;

					trueSpec = cvSpec;
					trueSpec.Scale("fullosc",0.0);
					std::cout<<"cv spec"<<std::endl;
					trueSpec.PrintCollapsedVector();
					appSpec = a_sinsqSpec[mi_base];
					// std::cout<<"app spec"<<std::endl;
					// appSpec.PrintCollapsedVector();
					appSpec.Scale("1e1p_fullosc",e_app);
					appSpec.Scale("bnb",-1.0*m_dis);
					appSpec.Scale("nue",-1.0*e_dis);
					appSpec.Scale("ccpi0",-1.0*m_dis);
					appSpec.Scale("ncpi0",-1.0*m_dis);
					appSpec.Scale("ext",0.0);
					trueSpec.Add(&appSpec);
					std::cout<<"oscillated spec"<<std::endl;
					trueSpec.PrintCollapsedVector();
					std::cout<<std::endl;
					if(mi_base==24 && uei_base==24 && umui_base ==24){
						for(int i = 0; i < nBins; i++){
							oscspecfile<<trueSpec.collapsed_vector[i]<<" ";
						}
						oscspecfile<<"\n";
					}

					// calculate the full covariance matrix based on the CV spectrum
					for(int i = 0; i < nBins; i++){
						for(int j = 0; j < nBins; j++){
								cov[i][j] = (*covFracSys_collapsed)[i][j]*cvSpec.collapsed_vector[i]*cvSpec.collapsed_vector[j];
								if(cov[i][j]==0){
									cov[i][j]=(*covFracSys_collapsed)[i][j];
								}

							if(i==j){
								// cov[i][j] += trueSpec.collapsed_vector[i];
								// add in stat errors start with CNP for "data" errors
								if (cvSpec.collapsed_vector[i] >0 ){
									cov[i][j] += 3.0 / (1.0/trueSpec.collapsed_vector[i] + 2.0/cvSpec.collapsed_vector[i]);
									// poisson version
									// cov[i][j] += pow(fakeData[i]-trueSpec.collapsed_vector[i],2) / (2.0*(trueSpec.collapsed_vector[i]- fakeData[i]+ fakeData[i]*log(fakeData[i]/trueSpec.collapsed_vector[i])));
								}
								else {
									cov[i][j] += cvSpec.collapsed_vector[i]/2.0;
								}
							}
						}//end of first bin loop
					}//end of second bin loop
					//save covar for plotting purposes
					if(mi_base==24 && uei_base==24 && umui_base ==24){
						for(short i = 0; i < nBins; i++){
							for(short j = 0; j < nBins; j++){
									covtotalfile<< cov[i][j]<<" ";
							}
							covtotalfile<<std::endl;
						}
					}
					// inv cov for chi2calc
					std::cout<<"determinanat: "<<cov.Determinant()<<std::endl;
					invcov = cov.Invert();
					if(mi_base==24 && uei_base==24 && umui_base ==24){
						for(short i = 0; i < nBins; i++){
							for(short j = 0; j < nBins; j++){
									covinvfile<< invcov[i][j]<<" ";
							}
							covinvfile<<std::endl;
						}
					}

					chisqTest = 0;
					for(int i = 0; i < nBins; i++){
						for(int j = 0; j < nBins; j++){
							// (obsi-predi)*(invcov)*(obsj-predj)
							chisqTest += (trueSpec.collapsed_vector[i] - cvSpec.collapsed_vector[i])*invcov[i][j]*(trueSpec.collapsed_vector[j] - cvSpec.collapsed_vector[j]);
						}
					}

					chifile << chisqTest<< std::endl;;

				}// end of loop over base umu
			}//end of loop over base ue
		}//end of loop over base mass
	}//end of not gen loop
	return 0;
}
