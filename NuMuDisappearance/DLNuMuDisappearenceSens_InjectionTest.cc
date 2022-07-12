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
	bool shapeonly = false;
	int miToDo = 0;
	bool doOneMass = false;
	int sin22ToDo = 0;
	bool doOneSin22 = false;
	bool statsOnly = false;
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
	const int dm2_grdpts(25), sin22th_grdpts(25);
	const int nFakeExp(50);
	const int nBins(19);
	double  mnu, sin22th, ue, um,_mnu, _sin22th, _ue, _um, chisqTest, chisqMin, chisqPT;
	int count;
	std::vector<double> fakeData;
	fakeData.resize(nBins);
	SBNspec truespec_first, appSpec_first, truespec_second, appSpec_second, testSpec;
	std::array<double,nBins>  a_pred;

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	if(!gen){
		// Load up the necesary bits and store them in vectors on the stack **
		std::string endStr = "_Bkg.SBNspec.root";
		SBNspec cvSpec(tag+endStr,xml);
		cvSpec.RemoveMCError();
		cvSpec.CollapseVector();


		std::vector<SBNspec> a_sinsqSpec;
		getPrecompSpec(a_sinsqSpec, dm2_grdpts, dm2_hibound, dm2_lowbound, xml, 1.0);



		// Bring in our covariance matrix!
		// Stats + sys
		// Load up cov matrix and add in detector variation component	**
		std::string covInFile;
		if (statsOnly) {covInFile = "JOSH_1m1p_StatOnly.SBNcovar.root";}
		else {covInFile = "JOSH_1m1p.SBNcovar.root";}
		TFile * fsys = new TFile(covInFile.c_str(),"read");
		TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
		TMatrixD * covFracSys_collapsed = (TMatrixD*)fsys->Get("frac_covariance_collapsed");

		std::cout<<"Loaded Covariance"<<std::endl;


		// // Load up oscillation model
		std::cout << "NUMU DISAPPEARANCE" << std::endl;

		float mnu_base, um_base, ue_base, sin22th_first;
		std::vector<std::vector<TMatrixD> > fracCovMats_grid = genCovMats(cvSpec,
																								a_sinsqSpec,
																								*covFracSys,
																								dm2_grdpts,
																								sin22th_grdpts,
																								sin22th_lowbound,
																								nBins
																							);

		std::vector<TH2D> gridPtloglike_h_v(nFakeExp,TH2D("ChiFit_Pts","ChiFit_Pts",25,0,25,25,0,25));



		// std::vector<double> injectionPtSin22th_v = {0.04,0.2,0.34,0.8};
		// double mVal = pow(2,0.5);
		// std::vector<double> injectionPtMNu_v     = {mVal,mVal,mVal,mVal};


		std::vector<double> injectionPtSin22th_v = {0.210};
		std::vector<double> injectionPtMNu_v     = {sqrt(3.374)};


		std::vector<double>    num_excluded_v(injectionPtMNu_v.size(),0);

		std::vector<std::vector<int>> bestFitsinPt_v(4,std::vector<int>(nFakeExp,0));
		std::vector<std::vector<int>> bestFitmiPt_v(4,std::vector<int>(nFakeExp,0));
		for (int inj_idx =0;inj_idx < injectionPtMNu_v.size(); inj_idx++){
			// if (inj_idx != 2){
			// 	continue;
			// }
			// Setting up a way to call the script piecemeal.
			//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
			// sin22th_first = pow(10.,((sin22thi_base+0.5)/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
			// mnu_base = pow(10.,((mi_base+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			sin22th_first = injectionPtSin22th_v[inj_idx];
			mnu_base      = injectionPtMNu_v[inj_idx];
			std::cout << "Running for Sin and Mass:" << sin22th_first << " " << mnu_base << "\n";
			// For Nue Appearance:
			// um_base = .135;	// global best fit
			// ue_base = pow(sin22th_first/float(4*um_base*um_base),.5);
			// For NuMu Disappearance:
			um_base = 1.0 ; //THIS VAUE IS ALSO DUMMY
			ue_base = 1.0 ; // Dummy Value to satisfy code
			// current test model
			std::cout << "NU Base Model: m41^2:" << mnu_base << " ue:" << ue_base << " sin2:" << sin22th_first << std::endl;

			// NuMu Disappearance Math described by equations 4.20-4.22 in Davio's thesis.
			// cvSpec.PrintCollapsedVector();
			// assert (1==2);

			truespec_first = cvSpec;
			bool doGen = true;
			SBNspec lostSpec = getAppSpec_exact(mnu_base,doGen);
			lostSpec.RemoveMCError();
			// continue;

			scaleDisappearedSpec(truespec_first, lostSpec, sin22th_first);
			std::cout<<"disappeared spec"<<std::endl;
			truespec_first.PrintCollapsedVector();
			std::cout << "Done\n";
			return 0;
			// makeSpectraPlots(cvSpec, appSpec_first, truespec_first);

			//  I'll be generating a new universe around the oscillated spectrum
			// std::cout<<"truespec size: "<<truespec_first.collapsed_vector.size()<<std::endl;
			// std::cout<<"matrix size:   "<<covFracSys->GetNrows()<<std::endl;

			 // I'll be generating a new universe around the oscillated spectrum
			SBNchi TrueChi(truespec_first,*covFracSys);
			// assert (1==2);
			TrueChi.ReloadCoreSpectrum(&truespec_first);
			TrueChi.pseudo_from_collapsed = true;
			TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomposition

			// Across several fake experiments:
			for(int expi = 0; expi < nFakeExp; expi++){
				std::cout << expi << " Experiment\n";
				std::vector<float> fakeData_f = TrueChi.GeneratePseudoExperiment();	// [fakedata] = vector of floats
				for (int bini = 0;bini<fakeData.size();bini++){
					fakeData[bini] = (double)fakeData_f[bini];
				}

				// print out new oscillated Spectrum
				std::cout<<"Oscillated Spectra At Injection Point:"<<std::endl;
				// truespec_first.PrintCollapsedVector();
				for(int i = 0; i < nBins; i++){
					std::cout<<truespec_first.collapsed_vector[i]<<" ,";
				}
				std::cout<<"\n";
				std::cout<<"Thrown Spectra: \n";
				for(int i = 0; i < nBins; i++){
					std::cout<<fakeData[i]<<" ,";
				}
				std::cout<<"\n";
				/*
				Surround with another grid search to find best fit
				original grid search variables:
				mi_base
				sin22thi_base
				*/
				double detTermTest, detTermPT, detTermMin;
				double chisqTest, chisqMin ,chisqPT;
				double halfRTermTest, halfRTermPT;
				double halfRTermMin = 99999999;
				double rTerm;

				int count=0;
				double sin22th_second;
				int min_sin22_idx=-1;
				int min_mi_idx=-1;
				double mnu_second;
				double best_mnu_second;
				double best_sin22th_second;
				for(int mi_second = 0; mi_second <dm2_grdpts; mi_second++){
					for(int sin22thi_second = 0; sin22thi_second <sin22th_grdpts; sin22thi_second++){
						TMatrixD cov(nBins,nBins), invcov(nBins,nBins), covCopy(nBins,nBins), isDiag(nBins,nBins);
						sin22th_second = pow(10.,((sin22thi_second+0.5)/float(sin22th_grdpts)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
						truespec_second = cvSpec;
						scaleDisappearedSpec(truespec_second, sbn::SBNspec(a_sinsqSpec[mi_second]), sin22th_second, false);


						if(shapeonly){
							scaleCovShapeOnly(cov, fracCovMats_grid[mi_second][sin22thi_second], fakeData, truespec_second.collapsed_vector, nBins); // cov, covfraccollapsed, obs, exp(gets modified), nbins
						}//End of shape only
						// Shape+Rate
						else{
							scaleCovShapeRate(cov, fracCovMats_grid[mi_second][sin22thi_second], fakeData, truespec_second.collapsed_vector, nBins); // cov, covfraccollapsed, obs, exp, nbins
						} //end of shape+rate
						// inv cov for chi2calc

						covCopy = cov;
						invcov = covCopy.Invert();

						calcRTerm(detTermTest, chisqTest, fakeData, truespec_second.collapsed_vector, invcov, cov, nBins); // obs, exp, inv, nbins
						halfRTermTest = detTermTest + chisqTest;
						if (inj_idx == 0){
							gridPtloglike_h_v[expi].SetBinContent(sin22thi_second+1,mi_second+1,halfRTermTest);
						}
						if ((mi_second == 0) && (0 == sin22thi_second)){
							detTermPT = detTermTest;
							chisqPT   = chisqTest;
							halfRTermPT = halfRTermTest;
						}
						if (halfRTermTest < halfRTermMin){
							chisqMin = chisqTest;
							detTermMin = detTermTest;
							halfRTermMin = halfRTermTest;
							min_sin22_idx = sin22thi_second;
							min_mi_idx = mi_second;
							best_mnu_second = mnu_second;
							best_sin22th_second = sin22th_second;
						}
					}//End second sin22 search
				}//End second m41 search
				std::cout << halfRTermPT << " " << halfRTermMin << " " << halfRTermPT - halfRTermMin << " HalfRTerms and diff\n";
				std::cout << "Min Sin22: "<< min_sin22_idx << " Min Mi:" << min_mi_idx << "\n";
				bestFitsinPt_v[inj_idx][expi] = min_sin22_idx;
				bestFitmiPt_v[inj_idx][expi]  = min_mi_idx;
				if ((halfRTermPT - halfRTermMin) > 5.90178){ // this is the crit rterm for the null pt.
					num_excluded_v[inj_idx]++;
				}
			}//End universe search
		}//End of Injection Points

		TFile outfile("injectionOutFile.root", "RECREATE");
		TTree gridTree("gridPtloglikes","gridPtloglikes");
		TH2D gridPtloglike_h("LogLikesFit_Pts_One","LogLikesFit_Pts_One",25,0,25,25,0,25);
		gridTree.Branch("gridPtloglike_h",&gridPtloglike_h);
		for (int idx =0; idx < gridPtloglike_h_v.size(); idx++){
			gridPtloglike_h = gridPtloglike_h_v[idx];
			gridTree.Fill();
		}
		gridTree.Write();
		outfile.Close();
		TGraph g("g","g");
		g.SetMarkerStyle(kFullStar);
		g.SetMarkerSize(5);
		g.SetMarkerColor(kBlack);
		for (int inj_idx =0 ; inj_idx<4;inj_idx++){
			gStyle->SetOptStat(0);
			std::string name = "BestFit_"+std::to_string(inj_idx);
			TCanvas can("canvas","canvas",1600,1200);
			TH2D bestFits(name.c_str(),name.c_str(),25,0,25,25,0,25);
			for (int expi=0; expi < nFakeExp;expi++){
				bestFits.Fill(bestFitsinPt_v[inj_idx][expi],bestFitmiPt_v[inj_idx][expi]);
			}
			bestFits.Draw("COLZ");
			std::vector<int> sidx_v = {7.5, 16.5, 19.5, 23.5};
			std::vector<int> midx_v = {14.5,14.5, 14.5, 14.5};
			g.SetPoint(0,sidx_v[inj_idx],midx_v[inj_idx]);
			g.Draw("PSAME");
			name = "BestFit_"+std::to_string(inj_idx)+".png";
			can.SaveAs(name.c_str());
		}

		std::cout <<"Ending Injection Point Test, Results:\n";
		for (int inj_idx = 0; inj_idx < injectionPtMNu_v.size(); inj_idx++){
			std::cout << "sin22: " << injectionPtSin22th_v[inj_idx] << " m2: " << injectionPtMNu_v[inj_idx] << " " << num_excluded_v[inj_idx] << " " << num_excluded_v[inj_idx]/nFakeExp << "\n";
		}
	}//End non-gen if statement
	std::cout << "\n";
	std::cout << "End of Script \n";
	return 0;
}
