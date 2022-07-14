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
		bkg.Scale("fullosc",0.0);
		bkg.WriteOut(tag+"_Bkg");

		float mnu;
		for(int mi = 0; mi < 100; mi++){
			mnu = pow(10.,(mi/float(100)*TMath::Log10(10./.1) + TMath::Log10(.1)));
			//std::cout << "MNU: " << mnu << std::endl;

			//Model: mnu, ue4, um4
			//we're precomputing, so we don't really care about the u's
			NeutrinoModel testModel(mnu, 1, 1);

			// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
			SBNgenerate * gen = new SBNgenerate(xml,testModel);
			// Write them to file
			gen->WritePrecomputedOscSpecs(tag);
		}
		return 1;
	}

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	if(!gen){
	
    // Create our unoscillated background
    SBNspec bkg(tag+"_Bkg.SBNspec.root",xml);
		
		//Bring in our covariance matrix!
		// Stats only.
		TMatrixD *cov;
		SBNchi uboone_chi_statsonly(bkg,true);

		// Stats + sys
		TFile * fsys = new TFile("DL.SBNcovar.root","read");
		cov = (TMatrixD*)fsys->Get("frac_covariance_DL");
		std::cout << "Frac Covariance Matrix" << std::endl;
	  	cov->Print();
		SBNchi uboone_chi(bkg,*cov);

		// Stats +  15% sys
/*		
		double flatSys = .15;
		TMatrixD flatCov;
		flatCov.ResizeTo(40,40);
		for(int i = 0; i < 40; i++){
			flatCov(i,i) = pow(flatSys,2);
		}
		flatCov.Print();
		SBNchi uboone_chi_flatsys(bkg,flatCov);
*/


		// Load up oscillation model
		NeutrinoModel nullModel(0,0,0);
		SBNosc osctrue(tag+"_Bkg.SBNspec.root",xml, nullModel);
		//SBNosc osctrue(tag+".SBNspec.root",xml, nullModel);
	
		// Print out a few sample spec
		//NeutrinoModel highdm2(sqrt(72.f),1.f,sqrt((1.f+sqrt(1.f-.2))/2.f));
		NeutrinoModel highdm2(sqrt(72.f),pow(.0015/4.f,.5),1.);
		NeutrinoModel meddm2(sqrt(2.f),pow(.006/4.f,.5),1.);
		
		SBNgenerate * gen1 = new SBNgenerate(xml,meddm2);
		SBNgenerate * gen2 = new SBNgenerate(xml,highdm2);
		gen1->WritePrecomputedOscSpecs(tag);
		gen2->WritePrecomputedOscSpecs(tag);

		SBNosc osc1 = osctrue;
		osc1.LoadModel(meddm2);
		osc1.OscillateThis(tag);		
		osc1.WriteOut("dm2_med");

		SBNosc osc2 = osctrue;
		osc2.LoadModel(highdm2);
		osc2.OscillateThis(tag);
		osc2.WriteOut("dm2_high");		

		if(numudis){
			// If we're doing numu disappearance:
			std::cout << "NUMU DISAPPEARANCE" <<  std::endl;
			float mnu, um, ue, sin22th;
			for(int mi = 0; mi < 100; mi++){
				for(int sin22thi = 0; sin22thi < 150; sin22thi++){

					sin22th = pow(10.,(sin22thi/float(100)*TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
					mnu = pow(10.,(mi/float(150)*TMath::Log10(10./.1) + TMath::Log10(.1)));
					ue = 1.f; // set sin22th(e e) = 0	
					um = sqrt((1.f+sqrt(1.f-sin22th))/2.f);
					NeutrinoModel testModel(mnu,ue,um);
				
					std::cout << "NU MODEL: " << mnu << " " << ue << " " << um << std::endl;
					SBNosc osc = osctrue;
					osc.SetBothMode();

					osc.LoadModel(testModel);
					osc.OscillateThis(tag);

					double chi2 = uboone_chi.CalcChi(&osc);
					double chi2_statsonly = uboone_chi_statsonly.CalcChi(&osc);
					//double chi2_flatsys = uboone_chi_flatsys.CalcChi(&osc);

					std::cout << "COUNT: " << (mi*100 + sin22thi)/float(100) << "%" << std::endl;
					std::cout << "ANS: " << pow(mnu,2) << " " << um << " " << 4*um*um*(1-um*um) << " " << chi2 << std::endl;
					std::cout << "ANS_STATSONLY: " << pow(mnu,2) << " " << um << " " << 4*um*um*(1-um*um) << " " << chi2_statsonly << std::endl;
					//std::cout << "ANS_FLATSYS: " << pow(mnu,2) << " " << um << " " << 4*um*um*(1-um*um) << " " << chi2_flatsys << std::endl;
				}
			}
		}
		else if(combined){
			std::cout << "COMBINED FIT" <<  std::endl;
			float mnu, um, ue;
			int mass_end;
			if(mass_start > 0){
				mass_end = mass_start + 1;
			}
			else{
				mass_start = 0;
				mass_end = 50;
			}

			for(int mi = mass_start; mi < mass_end; mi++){
				for(int umi = 0; umi < 25; umi++){
					for(int uei = 0; uei < 25; uei++){

						um = pow(10.,(umi/float(25)*TMath::Log10(1./1e-3) + TMath::Log10(1e-3)));
						ue = pow(10.,(uei/float(25)*TMath::Log10(1./1e-3) + TMath::Log10(1e-3)));
						mnu = pow(10.,(mi/float(50)*TMath::Log10(10./.1) + TMath::Log10(.1)));
						NeutrinoModel testModel(mnu,ue,um);
						std::cout << "NU MODEL: " << mnu << " " << ue << " " << um << std::endl;

						SBNosc osc = osctrue;
						osc.LoadModel(testModel);
						osc.OscillateThis(tag);

						if(umi == 7 && mi == 7 ){
							osc.WriteOut("osctest");
							bkg.WriteOut("bkgtest");
						}

						double chi2 = uboone_chi.CalcChi(&osc);
						double chi2_statsonly = uboone_chi_statsonly.CalcChi(&osc);
						std::cout << "COUNT: " << (mi*25*25 + umi*25 + uei)/float(50*25*25/100) << "%" << std::endl;
						std::cout << "ANS: " << pow(mnu,2) << " " << um << " " << ue << " " << chi2 << std::endl;
						std::cout << "ANS_STATSONLY: " << pow(mnu,2) << " " << um << " " << ue << " " << chi2_statsonly << std::endl;
					}
				}
			}
		}
		else{
			// If we're doing nue appearance
			std::cout << "NUE APPEARANCE" << std::endl;
			float mnu, um, ue, sin22th;
			for(int mi = 0; mi < 50; mi++){
				for(int sin22thi = 0; sin22thi < 50; sin22thi++){

					sin22th = pow(10.,(sin22thi/float(50)*TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
					mnu = pow(10.,(mi/float(50)*TMath::Log10(10./.1) + TMath::Log10(.1)));
					um = .09;//1.f;	// set sin22th(mu mu) = 0
					ue = pow(sin22th/float(4*um*um),.5);
					NeutrinoModel testModel(mnu,ue,um);

					std::cout << "NU MODEL: " << mnu << " " << ue << " " << um << std::endl;
					SBNosc osc = osctrue;

					osc.LoadModel(testModel);
					osc.OscillateThis(tag);

					double chi2 = uboone_chi.CalcChi(&osc);
					double chi2_statsonly = uboone_chi_statsonly.CalcChi(&osc);
					std::cout << "COUNT: " << (mi*50 + sin22thi)/float(50*.50) << "%" << std::endl;
					std::cout << "ANS: " << pow(mnu,2) <<  " " << ue << " " << um << " " << 4*um*um*ue*ue << " " << chi2 << std::endl;
					std::cout << "ANS_STATSONLY: " << pow(mnu,2) <<  " " << ue << " " << um << " " << 4*um*um*ue*ue << " " << chi2_statsonly << std::endl;
				}
			}
		}

		return 0;
	}
}
