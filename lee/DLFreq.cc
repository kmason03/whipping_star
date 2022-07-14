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
#include "SBNcovariance.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"
#include "prob.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

int main(int argc, char* argv[])
{
	std::string xml = "numu_disp.xml";
	int iarg = 0;
	opterr=1;
	int index;
	bool gen = false;
	bool numudis = false;
	bool combined = false;
	int mass_start = -1;

	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"gen",	no_argument, 0, 'g'},
		{"dis",	no_argument,0,'d'},
		{"comb", no_argument,0,'c'},
		{"part", required_argument,0,'p'},
		{0,			no_argument, 		0,  0},
	};

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:dscp:g", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'g':
				gen = true;
				break;
			case 'd':
				numudis = true;
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

	////////////////////////////////////////////////////////////////////////////////////////
	std::string tag = "DL";

	const double dm2_lowbound(0.01), dm2_hibound(100);
  const double sin22th_lowbound(1e-2), sin22th_hibound(1.0);
  const int dm2_grdpts(25), sin22th_grdpts(25);
  const int nFakeExp(1000);
	const int nBins(19);
	const bool shapeonly(false);

	std::string s_fracDetSysMatrix = "/uboone/app/users/dcianci/DLLEESensitivity/Thesis2020/whipping_star/data/detsys_frac_apr7.txt";
	////////////////////////////////////////////////////////////////////////////////////////
	double sin22th, mnu, ue, um, _mnu, _sin22th, _ue, _um, chisqTest, chisqMin, chisqPT;
	int count;
	TMatrixD cov(nBins,nBins), invcov(nBins,nBins);
	std::array<double,nBins> a_specTest, a_pred;
	std::vector<float> fakeData;
	fakeData.resize(nBins);
	SBNspec trueSpec, disSpec, testSpec;
	////////////////////////////////////////////////////////////////////////////////////////

	//PART 1: precompute a covariance matrix
	if(gen){

		SBNcovariance _covar(xml);
		_covar.FormCovarianceMatrix(tag);	
		_covar.PrintMatricies(tag);
		_covar.frac_covariance.Print();

		NeutrinoModel nullModel(0, 0, 0);
		SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
		SBNspec bkg = bkgo->spec_central_value;
		bkg.WriteOut(tag+"_Bkg");
		
		delete bkgo;

		for(int mi = 0; mi < dm2_grdpts; mi++){
			mnu = pow(10.,((mi+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			std::ostringstream out;
			NeutrinoModel testModel(mnu, 1, 1);
			SBNgenerate * gens = new SBNgenerate(xml,testModel);
			gens->WritePrecomputedOscSpecs(tag);
			delete gens;
		}
		return 1;
	}

	if(shapeonly){
		std::cout << "SHAPEONLY" << std::endl;
	}

	// Print binedges for easier plotting **
	std::cout << "BE_dm2: ";
	for(int mi = 0; mi <= dm2_grdpts; mi++){
		mnu = pow(10.,((mi)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		std::cout << pow(mnu,2) << " ";
	}
	std::cout << std::endl;

	std::cout << "BE_sin22th: ";
	for(int si = 0; si <= sin22th_grdpts; si++){					
		sin22th = pow(10.,((si)/float(sin22th_grdpts)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
		std::cout << sin22th << " ";
	}
	std::cout << std::endl;
	
	// Load up the necesary bits and store them in vectors on the stack **
  SBNspec cvSpec(tag+"_Bkg.SBNspec.root",xml);
	cvSpec.CollapseVector();
	std::array<SBNspec,dm2_grdpts> a_sinsqSpec;
	for(int mi = 0; mi < dm2_grdpts; mi++){
		mnu = pow(10.,((mi+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		std::stringstream stream;
		stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
		std::string infile = "DL_SINSQ_dm_"+stream.str()+".SBNspec.root";
		auto inspec = SBNspec(infile,xml);
		inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background. 
		inspec.CollapseVector();
		a_sinsqSpec[mi] = inspec;
	}

	// Load up cov matrix and add in detector variation component	**
	TFile * fsys = new TFile("DL.SBNcovar.root","read");
	TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
		
	std::ifstream file;
	double ddummy;
	file.open(s_fracDetSysMatrix);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << s_fracDetSysMatrix << std::endl;
	for(short i = 0; i < nBins; i++){
		for(short j = 0; j < nBins; j++){
			file >> ddummy;
			(*covFracSys)(i,j) += ddummy;
		}
	}
	file.close();

	SBNchi TrueChi(cvSpec,*covFracSys);

	// For each point in parameter space, pT (true point) **
	for(int pTdm2i = 0; pTdm2i < dm2_grdpts; pTdm2i++){
  	for(int pTsin22thi = 0; pTsin22thi < sin22th_grdpts; pTsin22thi++){

			sin22th = pow(10.,((pTsin22thi+0.5)/float(sin22th_grdpts)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
			mnu = pow(10.,((pTdm2i+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

			std::cout << "PT(sin22th,dm2): " << sin22th << " " << pow(mnu,2) << " --- " << pTsin22thi << " " << pTdm2i <<  std::endl;

		 	trueSpec = cvSpec;
			disSpec = a_sinsqSpec[pTdm2i];
			disSpec.Scale("bnb",-sin22th);
			trueSpec.Add(&disSpec);
			trueSpec.CollapseVector();

			TrueChi.ReloadCoreSpectrum(&trueSpec);
			TrueChi.pseudo_from_collapsed = true;
			TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomposition	
			
			// Across several fake experiments:
			std::cout << "PSEUDOEXP_DELTACHI2: ";
			for(int expi = 0; expi < nFakeExp; expi++){			
				fakeData = TrueChi.GeneratePseudoExperiment();	// [fakedata] = vector of floats

				// now find the best fit.
				chisqMin = 99999;	count=0;
				for(int mi = 0; mi < dm2_grdpts; mi++){
					for(int si = 0; si < sin22th_grdpts; si++){
							
						_sin22th = pow(10.,((si+0.5)/float(sin22th_grdpts)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
						for(int i = 0; i < nBins; i++)
							a_pred[i] = cvSpec.collapsed_vector[i] - _sin22th*a_sinsqSpec[mi].collapsed_vector[i];

						if(shapeonly){
							double covIntegral(0.0), predIntegral(0.0), obsIntegral(0.0),fnorm;
							for(int i = 0; i < nBins; i++){
								predIntegral+=a_pred[i];
								obsIntegral+=fakeData[i];
								for(int j = 0; j < nBins; j++){
									covIntegral += (*covFracSys)[i][j]*a_pred[i]*a_pred[j];
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
						else{
							for(int i = 0; i < nBins; i++){
								for(int j = 0; j < nBins; j++){
									cov[i][j] = (*covFracSys)[i][j]*a_pred[i]*a_pred[j];
									if(i==j)
										cov[i][j] += a_pred[i];
								}
							}
						}	
						invcov = cov.Invert();
					
						chisqTest = 0;
						for(int i = 0; i < nBins; i++){
							for(int j = 0; j < nBins; j++){
								chisqTest += (fakeData[i] - a_pred[i])*invcov[i][j]*(fakeData[j] - a_pred[j]);
							}
						}
						if(mi == pTdm2i && si == pTsin22thi)
							chisqPT = chisqTest;

						chisqMin = std::min(chisqTest,chisqMin);
						//count++;
					}
				}
				std::cout << chisqPT-chisqMin << " ";
			}
			std::cout << std::endl;
			return 1;
		}
	}
	return 1;
}
