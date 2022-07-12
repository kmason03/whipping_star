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
#include <assert.h>

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
#include "TDirectory.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"
#include "SBNcovariance.h"
#include "prob.h"
#include "MillsFunctions.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;


int main(int argc, char* argv[]){
	std::string xml = "/cluster/tufts/wongjiradlab/jmills09/whipping_star/xml/numudisappearance.xml";
	int iarg = 0;
	opterr=1;
	int index;
	bool gen = false;
	bool numudis = false;
	bool nueapp = true;
	bool combined = false;
	int mass_start = -1;
	bool shapeonly = true;
	int miToDo = 0;
	bool doOneMass = false;
	int sin22ToDo = 0;
	bool doOneSin22 = false;
	const struct option longopts[] ={
	{"xml", 		required_argument, 	0, 'x'},
	{"gen",	no_argument, 0, 'g'},
	{"dis",	no_argument,0,'d'},
	{"app", no_argument,0,'a'},
	{"comb", no_argument,0,'c'},
	{"part", required_argument,0,'p'},
	{"mibase", required_argument,0,'mi'},
	{"sinbase", required_argument,0,'si'},

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
			case 'si':
				doOneSin22 = true;
				sin22ToDo = atoi(optarg);
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
	const int nFakeExp(1000);
	const int nBins(19);
	const int dm2_grdpts(25), sin22th_grdpts(25);
	bool statsOnly = false;
	double  mnu, sin22th, ue, um,_mnu, _sin22th, _ue, _um, chisqTest, chisqMin, chisqPT;
	int count;
	std::vector<double> fakeData;
	fakeData.resize(nBins);
	SBNspec truespec_first, appSpec_first;

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	// save output to text files
	// Load up the necesary bits and store them in vectors on the stack **
	SBNspec cvSpec(tag+"_Bkg.SBNspec.root",xml);
	cvSpec.RemoveMCError();
	cvSpec.CollapseVector();
	// cvSpec.Scale("fullosc",0.0);
	// std::cout << "\n\n";
	// cvSpec.PrintCollapsedVector();
	// std::cout << "\n\n";
	// assert (1==2);
	// Setup Root file for output as well:


	std::vector<SBNspec> a_sinsqSpec;
	getPrecompSpec(a_sinsqSpec, dm2_grdpts, dm2_hibound, dm2_lowbound, xml, 1.0);
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

	std::vector<std::vector<TMatrixD> > fracCovMats_grid = genCovMats(cvSpec,
																							a_sinsqSpec,
																							*covFracSys,
																							dm2_grdpts,
																							sin22th_grdpts,
																							sin22th_lowbound,
																							nBins
																						);

	float mnu_base, um_base, ue_base, sin22th_first;
	// first get the central grid point we are testing

	int mi_base = 11; //Corresponds to value of 0.575
	int sin22thi_base = 23; // Corresponds to value of 0.69
	sin22th_first = pow(10.,((sin22thi_base+0.5)/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
	mnu_base = pow(10.,((mi_base+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

	std::cout << "Running for Sin and Mass:" << sin22th_first << " " << mnu_base << "\n";

	um_base = 1.0 ; //THIS VAUE IS ALSO DUMMY
	ue_base = 1.0 ; // Dummy Value to satisfy code
	// current test model
	std::cout << "NU Base Model: m41^2:" << pow(mnu_base,2) << " ue:" << ue_base << " sin2:" << sin22th_first << std::endl;
 	// NuMu Disappearance Math described by equations 4.20-4.22 in Davio's thesis.
	truespec_first = cvSpec;
	scaleDisappearedSpec(truespec_first, sbn::SBNspec(a_sinsqSpec[mi_base]), sin22th_first);
	std::cout<<"oscillated spec"<<std::endl;
	truespec_first.PrintCollapsedVector();

	//  I'll be generating a new universe around the oscillated spectrum
	SBNchi TrueChi(truespec_first,*covFracSys);
	TrueChi.ReloadCoreSpectrum(&truespec_first);
	TrueChi.pseudo_from_collapsed = true;
	TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomposition

	TH1F oscSpec("1000 Pseudo-Experiments","1000 Pseudo-Experiments",19,250,1200);
	oscSpec.SetXTitle("Reconstructed Neutrino Energy");
	oscSpec.SetYTitle("Number of Events");
	oscSpec.SetLineColor(2);
	oscSpec.SetLineWidth(6);
	for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
		oscSpec.SetBinContent(bin+1,truespec_first.collapsed_vector[bin]);
	}
	TMatrixD cov(nBins,nBins), invcov(nBins,nBins), covCopy(nBins,nBins);
	for (int i=0;i<nBins;i++){
		for (int j=0;j<nBins;j++){
			covCopy[i][j] = (*covFracSys_collapsed)[i][j]*truespec_first.collapsed_vector[i]*truespec_first.collapsed_vector[j];
		}
	}


	std::vector<double> sigmaVals(nBins,0);
	for (int i=0;i<nBins;i++){
		// Add sqrt(N) in quadrature to diagonal term
		sigmaVals[i] = sqrt((covCopy[i][i]*covCopy[i][i]));
		// sigmaVals[i] = sqrt((covCopy[i][i]*covCopy[i][i])+truespec_first.collapsed_vector[i]);
	}
	double max_bin_content;
	// One Sigma
	TH1F oneSigmaHigh("s1h","s1h",19,250,1200);
	oneSigmaHigh.SetLineColor(kBlack);
	oneSigmaHigh.SetLineWidth(6);
	oneSigmaHigh.SetLineStyle(kDashed);
	for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
		oneSigmaHigh.SetBinContent(bin+1,truespec_first.collapsed_vector[bin]+sqrt(sigmaVals[bin]));
	}
	TH1F oneSigmaLow("s1l","s1l",19,250,1200);
	oneSigmaLow.SetLineColor(kBlack);
	oneSigmaLow.SetLineWidth(6);
	oneSigmaLow.SetLineStyle(kDashed);
	for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
		oneSigmaLow.SetBinContent(bin+1,truespec_first.collapsed_vector[bin]-sqrt(sigmaVals[bin]));
	}
	// Two Sigma
	TH1F twoSigmaHigh("s2h","s2h",19,250,1200);
	twoSigmaHigh.SetLineColor(kBlack);
	twoSigmaHigh.SetLineWidth(6);
	twoSigmaHigh.SetLineStyle(9);
	for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
		twoSigmaHigh.SetBinContent(bin+1,truespec_first.collapsed_vector[bin]+(2*sqrt(sigmaVals[bin])));
	}
	TH1F twoSigmaLow("s2l","s2l",19,250,1200);
	twoSigmaLow.SetLineColor(kBlack);
	twoSigmaLow.SetLineWidth(6);
	twoSigmaLow.SetLineStyle(9);
	for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
		twoSigmaLow.SetBinContent(bin+1,truespec_first.collapsed_vector[bin]-(2*sqrt(sigmaVals[bin])));
	}
	// Three Sigma
	TH1F threeSigmaHigh("s3h","s3h",19,250,1200);
	threeSigmaHigh.SetLineColor(kBlack);
	threeSigmaHigh.SetLineWidth(6);
	threeSigmaHigh.SetLineStyle(kSolid);
	for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
		double val = truespec_first.collapsed_vector[bin]+(3*sqrt(sigmaVals[bin]));
		if (val > max_bin_content){
			max_bin_content = val;
		}
		threeSigmaHigh.SetBinContent(bin+1,val);
	}
	TH1F threeSigmaLow("s3l","s3l",19,250,1200);
	threeSigmaLow.SetLineColor(kBlack);
	threeSigmaLow.SetLineWidth(6);
	threeSigmaLow.SetLineStyle(kSolid);
	for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
		threeSigmaLow.SetBinContent(bin+1,truespec_first.collapsed_vector[bin]-(3*sqrt(sigmaVals[bin])));
	}

	std::vector<TH1F> fakeSpec_v(nFakeExp,TH1F("f0","f0",19,250,1200));
	// Across several fake experiments:
	std::vector<int> nOutside_v(19,0);
	int nOutside =0;
	for(int expi = 0; expi < nFakeExp; expi++){
		std::cout << expi << " Experiment\n";
		std::vector<float> fakeData_f = TrueChi.GeneratePseudoExperiment();	// [fakedata] = vector of floats
		for (int bini = 0;bini<fakeData.size();bini++){
			fakeData[bini] = (double)fakeData_f[bini];
		}
		// scaleCovShapeRate(covSWAP, fracCovMats_grid[mi_base][sin22thi_base], fakeData, truespec_first.collapsed_vector, nBins); // cov, covfraccollapsed, obs, exp, nbins

		std::cout<<"osc spec"<<std::endl;
		truespec_first.PrintCollapsedVector();
		std::cout<<"new uni spec: ";
		for(int i = 0; i < nBins; i++){
			std::cout<<fakeData[i]<<" ";
		}
		std::cout<<"\n";

		for (int bin=0; bin<fakeData.size(); bin++){
			if (fakeData[bin] > max_bin_content) max_bin_content = fakeData[bin];
			fakeSpec_v[expi].SetBinContent(bin+1,fakeData[bin]);
			if (fakeData[bin] > threeSigmaHigh.GetBinContent(bin+1)){
				nOutside++;
				nOutside_v[bin]++;
			}
			if (fakeData[bin] < threeSigmaLow.GetBinContent(bin+1)){
				nOutside++;
				nOutside_v[bin]++;
			}
		}
	}//End universe search



	gStyle->SetOptStat(0);
	std::vector<int> colorStarts = {kOrange, kPink, kMagenta, kViolet, kAzure, kTeal, kSpring};
	std::vector<int> colorTots;
	for (int i=0; i<colorStarts.size();i++){
		for (int j=0;j<20;j++){
			colorTots.push_back(colorStarts[i]+j-9);
		}
	}
	TLegend legend(0.7,0.70,0.90,0.9);
	legend.AddEntry(&oscSpec,"Disappeared CV Spectra");
	legend.AddEntry(&oneSigmaHigh,"1 Sigma");
	legend.AddEntry(&twoSigmaHigh,"2 Sigma");
	legend.AddEntry(&threeSigmaHigh,"3 Sigma");

	TCanvas can("c","c",1600,1000);
	oscSpec.Draw("");
	legend.Draw("SAME");
	oscSpec.SetMaximum(max_bin_content*1.05);
	for (int fi = 0;fi<nFakeExp;fi++){
		fakeSpec_v[fi].SetLineColor(colorTots[fi%colorTots.size()]);
		fakeSpec_v[fi].Draw("SAME");
	}
	oneSigmaHigh.Draw("SAME");
	twoSigmaHigh.Draw("SAME");
	threeSigmaHigh.Draw("SAME");
	oneSigmaLow.Draw("SAME");
	twoSigmaLow.Draw("SAME");
	threeSigmaLow.Draw("SAME");
	oscSpec.Draw("SAME");

	can.SaveAs("FakeDataDistribution.png");

	std::cout << "End of Script \n";
	return 0;
}
