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
  const int dm2_grdpts(25), sin22th_grdpts(15);
  const int nFakeExp(1000);
	const int nBins(19);

	std::string s_fracDetSysMatrix = "/uboone/app/users/dcianci/DLLEESensitivity/Thesis2020/whipping_star/data/detsys_frac_apr7.txt";
	////////////////////////////////////////////////////////////////////////////////////////
	double sin22th, mnu, ue, um, _mnu, _sin22th, _ue, _um, chisqTest, chisqMin, chisqPT;
	int count;
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
			mnu = pow(10.,((mi+0.5)/float(dm2_grdpts+1)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			std::ostringstream out;
			NeutrinoModel testModel(mnu, 1, 1);
			SBNgenerate * gens = new SBNgenerate(xml,testModel);
			gens->WritePrecomputedOscSpecs(tag);
			delete gens;
		}
		return 1;
	}


  SBNspec specCV(tag+"_Bkg.SBNspec.root",xml);

	// frontload the precomputing. normally we'd do this in a different stage, but at this scale it is much faster to do them all at once and store them on the stack


	std::vector<SBNspec> v_sinsqspec;
	v_sinsqspec.resize(dm2_grdpts);

	
	for(int mi = 0; mi < dm2_grdpts; mi++){
		mnu = pow(10.,((mi+0.5)/float(dm2_grdpts+1)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		std::stringstream stream;
		stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
		std::string infile = "DL_SINSQ_dm_"+stream.str()+".SBNspec.root";
		v_sinsqspec[mi] = SBNspec(infile,xml);
	}

	// Load up cov matrix and add in detector variation component	
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

	/*	
	for(short i = 0; i < 2*nBins; i++){			// we have low/no stats in some ext bins, so let's just manually patch that up for now
		for(short j = 0; j < 2*nBins; j++){
			if(i >= nBins || j >= nBins)
				(*covFracSys)(i,j)=0;
			if(i==j and i >= nBins)
				(*covFracSys)(i,j)=0.0001;
		}
	}
	*/
	

	std::array<double,19> a_specTest;
	SBNspec trueSpec, disSpec, testSpec;
	TMatrixD cov;

	// For each point in parameter space, pT (true point)
	for(int pTdm2i = 0; pTdm2i < dm2_grdpts; pTdm2i++){
  	for(int pTsin22thi = 0; pTsin22thi < sin22th_grdpts; pTsin22thi++){
  
			sin22th = pow(10.,((pTsin22thi+0.5)/float(sin22th_grdpts+1)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
			mnu = pow(10.,((pTdm2i+0.5)/float(dm2_grdpts+1)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

			std::cout << "PT(sin22th,dm2): " << sin22th << " " << pow(mnu,2) << " --- " << pTsin22thi << " " << pTdm2i <<  std::endl;

		 	trueSpec = specCV;
			disSpec = v_sinsqspec[pTdm2i];
			disSpec.ScaleAll(-sin22th);
			trueSpec.Add(&disSpec);
			
			SBNchi TrueChi(trueSpec,*covFracSys);
			TrueChi.pseudo_from_collapsed = true;
			TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomp. it's going to give an error but it's wrong...
	
		
			// Across several fake experiments:
			std::cout << "PSEUDOEXP_DELTACHI2: ";
			for(int expi = 0; expi < nFakeExp; expi++){			
				auto fakeData = TrueChi.GeneratePseudoExperiment();	// [fakedata] = vector of doubles

				// now find the best fit. and since we already precomputed across the grid, we'll just use that.
				chisqMin = 99999;	count=0;
				for(int mi = 0; mi < dm2_grdpts; mi++){
					for(int si = 0; si < sin22th_grdpts; si++){
						//std::cout << "Progress: " << float(count)/(dm2_grdpts*sin22th_grdpts)/(100.f)) << "\% \r";						
						
						_mnu = pow(10.,((mi+0.5)/float(dm2_grdpts+1)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
						_sin22th = pow(10.,((si+0.5)/float(sin22th_grdpts+1)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
						
						auto testSpec = specCV;
						auto testDisSpec = v_sinsqspec[pTdm2i];
						testDisSpec.ScaleAll(-_sin22th);
						testSpec.Add(&testDisSpec);

						SBNchi testChi(testSpec,*covFracSys);
						chisqTest = testChi.CalcChi(&fakeData);

						if(mi == pTdm2i && si == pTsin22thi)
							chisqPT = chisqTest;

						chisqMin = std::min(chisqTest,chisqMin);
						//count++;
					}
				}
				std::cout << chisqPT-chisqMin << " ";
			}
			std::cout << std::endl;
		}
	}
	return 1;
}
