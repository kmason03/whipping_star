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
#include "SBNcls.h"
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
	std::string xml = "example.xml";
	int iarg = 0;
	opterr=1;
	int index;
	bool sample_from_covariance = false;

	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"covariance", 		no_argument,0,'c'},
		{0,			no_argument, 		0,  0},
	};

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:c", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'c':
				sample_from_covariance = true;
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				std::cout<<"\t-c\t--covariance\t\tSample from covariance matrix instead of Poisson"<<std::endl;
				return 0;
		}
	}

	std::string tag = "EXAMPLE3";

    	SBNspec sig("EXAMPLE1.SBNspec.root",xml);
	sig.Scale("leesignal",1.5);
	
	SBNspec bkg("EXAMPLE1.SBNspec.root",xml);
	bkg.Scale("leesignal",0.0);

	// Stats + sys
	TFile * fsys = new TFile("EXAMPLE1.SBNcovar.root","read");
	TMatrixD * cov = (TMatrixD*)fsys->Get("frac_covariance_EXAMPLE1");

	SBNchi *chi = new SBNchi(bkg,*cov);
	SBNchi *chi_statonly = new SBNchi(bkg);

	SBNcls cls_factory(&sig, &bkg,*cov);
	if(sample_from_covariance) cls_factory.setSampleCovariance();

	//SBNspec test = chi->sampleCovariance(&sig);
	//test.writeOut("test.root");
	//return 0;

	int num_MC_events = 100000;
	cls_factory.calcCLS(num_MC_events, tag);

	return 0;
}
