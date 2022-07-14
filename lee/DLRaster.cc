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

	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"gen",	no_argument, 0, 'g'},
		{"dis",	no_argument,0,'d'},
		{"comb", no_argument,0,'c'},
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
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}
	}

	std::string tag = "DL";
	
	int nmnu(200), nsin22th(100);
	double mnu_low(.1), mnu_high(10.);
	double sin22th_low(1e-5), sin22th_high(1.0);

	//PART 1: precompute all of our sin and sin2 amplitudes so we don't need to later
	if(gen){

		NeutrinoModel nullModel(0, 0, 0);
		SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
		SBNspec bkg = bkgo->spec_central_value;
		bkg.WriteOut(tag+"_Bkg");

		float mnu;
		for(int mi = 0; mi < nmnu; mi++){
			mnu = pow(10.,(mi/double(nmnu)*TMath::Log10(mnu_high/mnu_low) + TMath::Log10(mnu_low)));

			//Model: mnu, ue4, um4
			//we're precomputing, so we don't really care about the u's
			NeutrinoModel testModel(mnu, .1, .1);

			// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
			SBNgenerate * gen = new SBNgenerate(xml,testModel);
			// Write them to file
			gen->WritePrecomputedOscSpecs(tag);

			delete gen;
		}
		return 1;
	}

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	
  // Create our unoscillated background
  SBNspec bkg(tag+"_Bkg.SBNspec.root",xml);
		
	//Bring in our covariance matrix!
	// Stats only.
	TMatrixD *cov;
	SBNchi uboone_chi_statsonly(bkg,true);

	// Stats + sys
	TFile * fsys = new TFile("DL.SBNcovar.root","read");
	cov = (TMatrixD*)fsys->Get("frac_covariance");
	std::cout << "Frac Covariance Matrix" << std::endl;
	 	cov->Print();
	SBNchi uboone_chi(bkg,*cov);

	// Load up oscillation model
	NeutrinoModel nullModel(0,0,0);
	SBNosc osctrue(tag+"_Bkg.SBNspec.root",xml, nullModel);

	std::cout << "NUMU DISAPPEARANCE" <<  std::endl;
	double mnu, um, ue, sin22th, chi2min, chi2min_statsonly, rasterpoint, rasterpoint_statsonly;
	for(int mi = 0; mi < nmnu; mi++){

		mnu = pow(10.,(mi/double(nmnu)*TMath::Log10(mnu_high/mnu_low) + TMath::Log10(mnu_low)));
		std::cout << "i, dm2: " << mi << ", " << pow(mnu,2) << std::endl;
		// find chi2min for raster scanning
		chi2min = 999999; chi2min_statsonly = 99999;
		rasterpoint = -1; rasterpoint_statsonly = -1;
		for(int sin22thi = 0; sin22thi < nsin22th; sin22thi++){

			sin22th = pow(10.,(sin22thi/double(nsin22th)*TMath::Log10(sin22th_high/sin22th_low) + TMath::Log10(sin22th_low)));
			ue = .001;	// this is arbitrary	
			um = sqrt((1.f+sqrt(1.f-sin22th))/2.f);
			NeutrinoModel testModel(mnu,ue,um);
				
			std::cout << "NU MODEL: " << mnu << " " << ue << " " << um << std::endl;
			SBNosc osc = osctrue;
			osc.SetBothMode();

			osc.LoadModel(testModel);
			osc.OscillateThis(tag);

			double chi2 = uboone_chi.CalcChi(&osc);
			double chi2_statsonly = uboone_chi_statsonly.CalcChi(&osc);

			std::cout << "COUNT: " << (mi*nmnu + sin22thi)/(double(nmnu*nsin22th)/100.0) << "%" << std::endl;
			std::cout << "ANS: " << pow(mnu,2) << " " << um << " " << sin22th << " " << chi2 << std::endl;
			std::cout << "ANS_STATSONLY: " << pow(mnu,2) << " " << um << " " << sin22th << " " << chi2_statsonly << std::endl;

			if(chi2 < chi2min)
				chi2min = chi2;
			else if(rasterpoint < 0 && chi2 - chi2min > 2.71)
				rasterpoint = sin22th;
					
			if(chi2_statsonly < chi2min_statsonly)
				chi2min_statsonly = chi2_statsonly;
			else if(rasterpoint_statsonly < 0 && chi2_statsonly - chi2min_statsonly > 2.71)
				rasterpoint_statsonly = sin22th;

			if((rasterpoint_statsonly > 0 && rasterpoint > 0) || sin22thi == nsin22th-1){
				if(rasterpoint_statsonly < 0)	rasterpoint_statsonly = sin22th_high;
				if(rasterpoint < 0) rasterpoint = sin22th_high;
				std::cout << "RASTER: " << pow(mnu,2) << " " << rasterpoint << std::endl;
				std::cout << "RASTER_STATSONLY: " << pow(mnu,2) << " " << rasterpoint_statsonly << std::endl;
				break;
			}
		}
	}

	return 0;
}
