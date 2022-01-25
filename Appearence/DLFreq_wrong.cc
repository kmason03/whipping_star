// In this version, I incorrectly build the covariance matrices for comparisons. it is kept for posterity until the new version is totally fixed.

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

	std::string tag = "DL";

	//PART 1: precompute all of our sin and sin2 amplitudes so we don't need to later
	if(gen){

		NeutrinoModel nullModel(0, 0, 0);
		SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
		SBNspec bkg = bkgo->spec_central_value;
		//SBNspec bkg("DL.SBNspec.root",xml);
		bkg.WriteOut(tag+"_NoOsc");
		
		delete bkgo;

		float mnu;
		for(int mi = 0; mi < 100; mi++){
			mnu = pow(10.,((mi+0.5)/float(101)*TMath::Log10(10./.1) + TMath::Log10(.1)));
			//std::cout << "MNU: " << mnu << std::endl;

			//Model: mnu, ue4, um4
			//we're precomputing, so we don't really care about the u's
			NeutrinoModel testModel(mnu, 1, 1);

			// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
			SBNgenerate * gens = new SBNgenerate(xml,testModel);
			// Write them to file
			gens->WritePrecomputedOscSpecs(tag);
			delete gens;
		}
		return 1;
	}

	//PART  2: Now that sin and sin2 libs are generated, do the rest...
    
	TRandom3 RanGen();

	double alpha(0.9);

	double dm2_lowbound(0.01), dm2_hibound(100);
  double sin22th_lowbound(1e-4), sin22th_hibound(1.0);
  int dm2_grdpts(2), sin22th_grdpts(2);

	int bfCourseness(10);
  int nFakeExp(50);
	int nChains(20),chainDepth(100),chainCt;
	double step,temp,chi2min;

	std::vector<double> v_chi2min;
	v_chi2min.resize(nFakeExp);

	// Load no osc spectrum
  SBNspec NoOsc("DL_NoOsc.SBNspec.root",xml);
	
	// Load up cov matrix	
	TFile * fsys = new TFile("DL.SBNcovar.root","read");
	auto cov = (TMatrixD*)fsys->Get("frac_covariance");
	std::cout << "Frac Covariance Matrix" << std::endl;	


	// Here, we can add in the detector systematic component to the covariance matrix
	//
	//



	double sin22th, mnu, ue, um, _mnu, _sin22th, _ue, _um, testChi2;
	int count;
	SBNosc TrueOsc(NoOsc);

	
	
	// Create an osc model for the no osc hypothesis. This will correspond to the "true" signal

	std::cout << "Let's GOOOOOOOOOO" << std::endl;
  // For each point in parameter space, pT (true point)
	for(int pTdm2i = 0; pTdm2i < dm2_grdpts; pTdm2i++){
  	for(int pTsin22thi = 0; pTsin22thi < sin22th_grdpts; pTsin22thi++){
  
			sin22th = pow(10.,((pTsin22thi+0.5)/float(sin22th_grdpts+1)*TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
			mnu = pow(10.,((pTdm2i+0.5)/float(dm2_grdpts+1)*TMath::Log10(10./.1) + TMath::Log10(.1)));
			ue = 0.f; // set sin22th(e e) = 0	
			um = sqrt((1.f+sqrt(1.f-sin22th))/2.f);
			NeutrinoModel trueModel(mnu,ue,um);

			TrueOsc.LoadModel(trueModel);
			TrueOsc.OscillateThis(tag);

			SBNchi TrueChi(TrueOsc,*cov)
			auto trueCov = TrueChi.m_matrix_systematics_collapsed;
			// get that collapsed cov matrix out
			//			gives tmatrix			GetCollapsed			 (note, this contains stat error, take that out) (re-add prediction stat error)
			//true_chi.m_matrix_systematics_collapsed   ( no stats error ) <--- use this one. easier.
			// Across several fake experiments:
      for(int expi = 0; expi < nFakeExp; expi++){
				std::cout << "Fake Data " << expi << std::endl;
				
				auto fakeData = TrueChi.GeneratePseudoExperiment();

				// now find the best fit. and since we already precomputed across the grid, we'll just use that.
				chi2min = 99999;	count=0;
				v_chi2min[expi] = chi2min;	// set to garbage number to start
				for(int mi = 0; mi < bfCourseness; mi++){
					for(int si = 0; si < bfCourseness; si++){
						std::cout << "Progress: " << float(count)/(pow(bfCourseness,2)/(100.f)) << "\% \r";						
		
						_mnu = pow(10.,(mi/float(bfCourseness)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
						_sin22th = pow(10.,(si/float(bfCourseness)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
						_ue = 1.f; // set sin22th(e e) = 0	
						_um = sqrt((1.f+sqrt(1.f-sin22th))/2.f);
						NeutrinoModel testModel(_mnu,_ue,_um);

						SBNosc osc(NoOsc);
						osc.LoadModel(testModel);
						auto pred = osc.Oscillate(tag);
						
						SBNchi PredChi(pred,*cov);
						auto predCov = PredChi.m_matrix_systematics_collapsed;

						// Add stat error for prediction;
						for(int i = 0; i < pred.size(); i++)		
							predCov(i,i) += pred[i];
						auto invPredCov = predCov.Invert();

						// Calculate Chi2
						testChi2 = 0;
						for(int i =0; i < pred.size();i++)
							for(int j = 0; j < pred.size();j++)
								testChi2 += (pred[i]-fakeData[i])*invPredCov(i,j)*(pred[j]-fakeData[j]);

						//std::cout << "tchi: " << testChi2 << std::endl;
						chi2min = std::min(testChi2,chi2min);
						count++;
					}
				}
				v_chi2min[expi] = chi2min;
				std::cout << std::endl << "chi2min " << expi << " = " << chi2min << std::endl;
			}

			// Now find effective chisq for... let's do 90%. We'll do this stupidly.
			double chi2_eff;
			for(int i = 0; i < nFakeExp; i++){
				// want to find the xth-lowest chi2min;
				auto _mit = std::min_element(v_chi2min.begin(),v_chi2min.end());
				auto _mindex = std::distance(v_chi2min.begin(), _mit);
				if(i/float(nFakeExp) < (1-alpha))
					v_chi2min.at(_mindex) = 9999;
				else{
					chi2_eff = v_chi2min.at(_mindex);
					break;
				}
			}
			std::cout << "Chi2Eff (" << pow(mnu,2) << ", " << sin22th << ") = " << chi2_eff << std::endl;
		}
	}

	return 1;









/*




		
	//Bring in our covariance matrix!
	// Stats only.
	TMatrixD *cov;
	SBNchi uboone_chi_statsonly(bkg,true);


	

	SBNchi uboone_chi(bkg,*cov);

	// Load up oscillation model
	NeutrinoModel nullModel(0,0,0);
	SBNosc osctrue(tag+"_Bkg.SBNspec.root",xml, nullModel);

	/*
	// Print out a few sample spec
	//NeutrinoModel highdm2(sqrt(72.f),1.f,sqrt((1.f+sqrt(1.f-.2))/2.f));
	NeutrinoModel highdm2(sqrt(72.f),1,.3);
		
	SBNgenerate * gen2 = new SBNgenerate(xml,highdm2);
	gen2->WritePrecomputedOscSpecs(tag);

	
	SBNosc osc2 = osctrue;
	osc2.SetBothMode();
	osc2.LoadModel(highdm2);
	osc2.OscillateThis(tag);
	osc2.WriteOut("dm2_high");		

	std::cout << "hidmchi: " << uboone_chi.CalcChi(&osc2) << std::endl;

	

	// If we're doing numu disappearance:
	std::cout << "NUMU DISAPPEARANCE" <<  std::endl;
	float mnu, um, ue, sin22th;
	for(int mi = 0; mi < 100; mi++){
		for(int sin22thi = 0; sin22thi < 100; sin22thi++){

				
			std::cout << "NU MODEL: " << mnu << " " << ue << " " << um << std::endl;
			SBNosc osc = osctrue;
			osc.SetBothMode();
			osc.LoadModel(testModel);
			osc.OscillateThis(tag);

			double chi2 = uboone_chi.CalcChi(&osc);
			double chi2_statsonly = uboone_chi_statsonly.CalcChi(&osc);
			//double chi2_flatsys = uboone_chi_flatsys.CalcChi(&osc);

			std::cout << "COUNT: " << (mi*100 + sin22thi)/float(100) << "%" << std::endl;
			std::cout << "ANS: " << pow(mnu,2) << " " << um << " " << sin22th << " " << chi2 << std::endl;
			std::cout << "ANS_STATSONLY: " << pow(mnu,2) << " " << um << " " << sin22th << " " << chi2_statsonly << std::endl;
			//std::cout << "ANS_FLATSYS: " << pow(mnu,2) << " " << um << " " << 4*um*um*(1-um*um) << " " << chi2_flatsys << std::endl;
		}
	}
	*/
}
