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

// void makeSpectraPlotsComp_Data(const std::vector<double>& cv, const std::vector<double>& disap, const std::vector<double>& best, int set, const std::vector<double>& unc);
// void makeSpectraPlotsComp_Data(const std::vector<double>& cv, const std::vector<double>& disap, const std::vector<double>& best, int set, const std::vector<double>& unc){
// 	gStyle->SetOptStat(0);
//   int nBins = 19;
// 	// std::string titleStr = "#nu_{mu} Dis. Fake Dataset: "+std::to_string(set);
//   std::string titleStr = "Muon Neutrino Disappearance Fake Dataset: "+std::to_string(set);
//   if (set == 3){
//     titleStr = "DL Muon Neutrino Disappearance, Data: 6.67e20";
//   }
//   TH1F label_h(titleStr.c_str(),titleStr.c_str(),19,250,1200);
// 	TH1F cvSpec_h("tmp0","tmp0",19,250,1200);
// 	TH1F disSpec_h("tmp","tmp",19,250,1200);
//   TH1F best_h("tmp2","tmp2",19,250,1200);
//
// 	cvSpec_h.SetLineColor(kBlue);
// 	disSpec_h.SetLineColor(kBlack);
//   best_h.SetLineColor(kRed);
//   cvSpec_h.SetLineWidth(5);
//   disSpec_h.SetLineWidth(5);
//   best_h.SetLineWidth(5);
//   disSpec_h.SetMarkerStyle(20);
//
//
//   TGraphErrors errors_cv;
//   errors_cv.SetLineColor(kBlue-9);
//   errors_cv.SetLineWidth(2.5);
//   errors_cv.SetFillStyle(3353);
//   gStyle->SetHatchesLineWidth(2.5);
//
//   errors_cv.SetFillColor(kBlue-9);
//
// 	// disSpec_h.SetLineStyle(kDashed);
// 	double maxVal = 0;
// 	double minVal = 0;
// 	for (int i = 0; i<19; i++){
// 		if (cv[i] > maxVal) maxVal = cv[i];
// 		// if (disap[i] > maxVal) maxVal = disap[i];
// 		if (cv[i] < minVal) minVal = cv[i];
// 		// if (disap[i] < minVal) minVal = disap[i];
// 		cvSpec_h.SetBinContent(i+1,cv[i]);
// 		disSpec_h.SetBinContent(i+1,disap[i]);
//     best_h.SetBinContent(i+1,best[i]);
//
//     disSpec_h.SetBinError(i+1,pow(disap[i],0.5));
//     errors_cv.SetPoint(i,275+i*50,cv[i]);
//     errors_cv.SetPointError(i,25,unc[i]);
//
//     // cvSpec_h.SetBinError(i+1, unc[i]);
// 	}
//   TCanvas c1("c","c",1200,1000);
//   double botfrac = 0.23;
//   double topfrac = botfrac + 0.01;
//   double rightmargin = 0.12;
//   double leftmargin = 0.12;
//
//   TPad  pad("upper_pad","",0,topfrac,1.,1.);
//   TPad rpad("ratio_pad","",0,0,1.,botfrac);
//   pad.SetFrameLineWidth(2);
//   pad.SetTopMargin(0.88*1/topfrac);
//   pad.SetBottomMargin(0.0);
//   pad.SetLeftMargin(leftmargin);
//   pad.SetRightMargin(rightmargin);
//   rpad.SetFrameLineWidth(2);
//   rpad.SetTopMargin(0.05);
//   rpad.SetBottomMargin(0.12*1/botfrac);
//   rpad.SetLeftMargin(leftmargin);
//   rpad.SetRightMargin(rightmargin);
//
//   pad.Draw();
//   pad.cd();
//   // TH1F frame = pad.DrawFrame(xmin,ymin,xmax,ymax)
//   c1.cd();
//   pad.cd();
//
//   label_h.SetXTitle("Neutrino Energy");
//   label_h.SetYTitle("Events");
// 	label_h.SetMaximum(maxVal*1.2);
//   disSpec_h.SetMaximum(maxVal*1.2);
//   label_h.Draw();
//   errors_cv.Draw("SAME 2");
// 	cvSpec_h.Draw("SAME");
//
// 	disSpec_h.Draw("SAMEP");
//   best_h.Draw("SAME");
//
//   TLegend legend(0.45,0.75,1-rightmargin,0.9);
// 	std::string cv_str  = "Null Model Spectra";
// 	std::string dis_str = "Fake Dataset "+std::to_string(set);
//   std::string best_sr = "Best Fit Disappearance Spectra";
//   if (set ==3){
//     cv_str  = "Null Oscillation Model";
//     dis_str = "6.67e20 Data";
//     best_sr = "Best Fit 3+1 Disappearance Model";
//   }
//   legend.AddEntry(&cvSpec_h, cv_str.c_str());
// 	legend.AddEntry(&disSpec_h,dis_str.c_str());
//   legend.AddEntry(&best_h,best_sr.c_str());
// 	legend.Draw();
//
//   c1.cd();
//   rpad.SetBorderSize(1);
//   rpad.Draw();
//   rpad.cd();
//   double minrat =0.0;
//   double maxrat =1.2;
//   if (set == 3){
//     minrat = 0.6;
//     maxrat = 1.4;
//   }
//   TH1F rframe = (*rpad.DrawFrame(250,minrat,1200,maxrat));
//
//   rframe.GetYaxis()->SetTitle("Ratio to Null");
//   rframe.GetXaxis()->SetLabelFont(43);
//   rframe.GetXaxis()->SetLabelSize(25);
//   rframe.GetYaxis()->SetLabelFont(43);
//   rframe.GetYaxis()->SetLabelSize(25);
//   rframe.GetYaxis()->SetNdivisions(7);
//   if (set == 3){
//     rframe.GetYaxis()->SetNdivisions(5);
//   }
//   rframe.GetXaxis()->SetTitleFont(43);
//   rframe.GetXaxis()->SetTitleSize(25);
//   rframe.GetYaxis()->SetTitleFont(43);
//   rframe.GetYaxis()->SetTitleSize(25);
//   rframe.GetYaxis()->SetTitleOffset(1.5);
//   // rframe.GetXaxis()->SetTitle("Q2 Reco");
//   rframe.GetXaxis()->SetTitle("Neutrino Energy");
//   rframe.GetXaxis()->SetTitleOffset(5);
//   rframe.Draw();
//   // TH1F cvSpec_ratio("tmp3","tmp3",nBins,0,60000);
//   // TH1F disSpec_ratio("tmp4","tmp4",nBins,0,60000);
//   // TH1F fakedata_ratio("tmp5","tmp5",nBins,0,60000);
//   TH1F label_ratio("tmp1","tmp1",nBins,250,1200);
//   TH1F cvSpec_ratio("tmp3","tmp3",nBins,250,1200);
//   TH1F best_ratio("tmp4","tmp4",nBins,250,1200);
//   TH1F disap_ratio("tmp5","tmp5",nBins,250,1200);
//   cvSpec_ratio.SetLineColor(kBlue);
//   best_ratio.SetLineColor(kRed);
//   disap_ratio.SetLineColor(kBlack);
//   cvSpec_ratio.SetLineWidth(5);
//   best_ratio.SetLineWidth(5);
//   disap_ratio.SetLineWidth(5);
//   disap_ratio.SetMarkerStyle(20);
//
//
//   TGraphErrors errors_rat_cv;
//   errors_rat_cv.SetLineColor(kBlue-9);
//   errors_rat_cv.SetLineWidth(2.5);
//   errors_rat_cv.SetFillStyle(3353);
//   errors_rat_cv.SetFillColor(kBlue-9);
//
//   for (int i=0;i<nBins;i++){
//     cvSpec_ratio.SetBinContent(i+1,1.0);
//     best_ratio.SetBinContent(i+1,best_h.GetBinContent(i+1)*1.0/cvSpec_h.GetBinContent(i+1));
//     disap_ratio.SetBinContent(i+1,disSpec_h.GetBinContent(i+1)*1.0/cvSpec_h.GetBinContent(i+1));
//
//     errors_rat_cv.SetPoint(i,275+i*50,1.0);
//     errors_rat_cv.SetPointError(i,25,unc[i]/cvSpec_h.GetBinContent(i+1));
//     disap_ratio.SetBinError(i+1,pow(disap[i],0.5)/disap[i]);
//
//     // cvSpec_ratio.SetBinError(i+1, unc[i]/cvSpec_h.GetBinContent(i+1));
//   }
//   label_ratio.Draw("SAME");
//   errors_rat_cv.Draw("SAME 2");
//   cvSpec_ratio.Draw("SAME");
//   best_ratio.Draw("SAME");
//   disap_ratio.Draw("SAMEP");
//
// 	std::string savestr = "data_spectras/spec_d"+std::to_string(set)+".png";
// 	c1.SaveAs(savestr.c_str());
// }



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
	const int nFakeExp(3);
	const int nBins(19);
	bool statsOnly = false;
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

		cvSpec.PrintFullVector();
		std::cout << "Here\n";
		// for (int abc = 0;abc<19;abc++){
		// 	std::cout << cvSpec.collapsed_vector[abc] << "\n";
		// }
		// assert (1==2);

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
		// first get the central grid point we are testing


		std::vector<TH2D> gridPtloglike_h_v(nFakeExp,TH2D("loglikeFit_Pts","loglikeFit_Pts",25,0,25,25,0,25));

		std::vector<std::vector<int>> bestFitsinPt_v(1,std::vector<int>(nFakeExp,0));
		std::vector<std::vector<int>> bestFitmiPt_v(1,std::vector<int>(nFakeExp,0));

		// Across several fake experiments:
		for(int expi = 0; expi < nFakeExp; expi++){
			// if (expi != 2){
			// 	continue;
			// }
			std::cout << expi << " Experiment\n";
			if (expi == 0){
				fakeData = {18,128,201,279,236,312,328,314,229,194,177,155,144,148,121,89,63,37,48};//fset1 scaled to datapot, and rounded
			}
			else if (expi == 1){
				fakeData = {23,185,292,335,353,400,462,435,372,305,279,229,190,190,134,119,75,60,54};//fset2 scaled to datapot, and rounded
			}
			else if (expi == 2){
				// THIS IS THE REAL DATA SELECTION DO NOT USE
				fakeData = {26.0, 192.0, 276.0 ,401.0, 389.0, 463.0, 439.0, 482.0, 395.0, 353.0, 303.0, 233.0, 240.0, 177.0, 118.0, 109.0, 109.0, 85.0, 58.0};
			}

			std::cout<<"\n";
			std::cout<<"'Data' Spectra: \n";
			for(int i = 0; i < nBins; i++){
				std::cout<<fakeData[i]<<" ,";
			}
			std::cout<<"\n";
			// mcstat error - found in selection grid from python, in PlottingScripts.py 1/np.sqrt(yerrsq_mc) for each bin)
			// std::vector<float> stkerr={0.4, 0.39555943, 0.14108268, 0.12801097, 0.10042527, 0.11516246, 0.1107993, 0.10101455, 0.09678998, 0.11705153, 0.11300651, 0.11581347, 0.14006509, 0.13711131, 0.1510809, 0.17920406, 0.17735261, 0.24390278, 0.26385786, 0.27852506};
			std::vector<float> stkerr={0.39555943, 0.14108268, 0.12801097, 0.10042527, 0.11516246, 0.1107993, 0.10101455, 0.09678998, 0.11705153, 0.11300651, 0.11581347, 0.14006509, 0.13711131, 0.1510809, 0.17920406, 0.17735261, 0.24390278, 0.26385786, 0.27852506};

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
			std::vector<double> bestFitSpec_v(nBins,0);
			for(int mi_second = 0; mi_second <dm2_grdpts; mi_second++){
				for(int sin22thi_second = 0; sin22thi_second <sin22th_grdpts; sin22thi_second++){
					// std::cout << mi_second << " " << sin22thi_second <<"\n";
					TMatrixD cov(nBins,nBins), invcov(nBins,nBins), covCopy(nBins,nBins), isDiag(nBins,nBins);
					sin22th_second = pow(10.,((sin22thi_second+0.5)/float(sin22th_grdpts)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
					mnu_second      = pow(10.,((mi_second+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

					if (mi_second ==0){
						std::cout << sin22thi_second << " " << sin22th_second << "\n";
					}
					truespec_second = cvSpec;
					scaleDisappearedSpec(truespec_second, sbn::SBNspec(a_sinsqSpec[mi_second]), sin22th_second, false);
					// truespec_second.PrintCollapsedVector();
					if(shapeonly){
						scaleCovShapeOnly(cov, fracCovMats_grid[mi_second][sin22thi_second], fakeData, truespec_second.collapsed_vector, nBins); // cov, covfraccollapsed, obs, exp(gets modified), nbins
					}//End of shape only
					// Shape+Rate
					else{
						scaleCovShapeRate(cov, fracCovMats_grid[mi_second][sin22thi_second], fakeData, truespec_second.collapsed_vector, nBins); // cov, covfraccollapsed, obs, exp, nbins
					} //end of shape+rate

					covCopy = cov;
					invcov = covCopy.Invert();
					// isDiag.Mult(invcov,cov);

					calcRTerm(detTermTest, chisqTest, fakeData, truespec_second.collapsed_vector, invcov, cov, nBins); // obs, exp, inv, nbins
					halfRTermTest = detTermTest + chisqTest;
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
						for (int iii=0; iii<nBins;iii++){
							bestFitSpec_v[iii] = truespec_second.collapsed_vector[iii];
						}

					}
					gridPtloglike_h_v[expi].SetBinContent(sin22thi_second+1,mi_second+1,halfRTermTest);
				}//End second sin22 search
			}//End second m41 search

			std::vector<double> unc_v(nBins,0);
			for (int i=0; i< nBins; i++){
				unc_v[i] = pow(fracCovMats_grid[0][0][i][i]*cvSpec.collapsed_vector[i]*cvSpec.collapsed_vector[i],0.5);
			}

			makeSpectraPlotsComp_Data(cvSpec.collapsed_vector, fakeData, bestFitSpec_v, expi+1, unc_v);

			std::cout << detTermPT << " " << chisqPT<< " "  << detTermMin << " " << chisqMin << " " << detTermPT - detTermMin << " " << chisqPT - chisqMin << " " << detTermPT - detTermMin + chisqPT - chisqMin  << " Chisquares\n";
			std::cout << "Min Sin22: "<< min_sin22_idx << " Min Mi:" << min_mi_idx << "\n";
			std::cout << "Best Sin22: "<< best_sin22th_second << " Best Mi:" << best_mnu_second << "\n";


			bestFitsinPt_v[0][expi] = min_sin22_idx;
			bestFitmiPt_v[0][expi]  = min_mi_idx;
			std::string histName = "uni_"+ZeroPadNumber(expi,4)+"_dChi2_"+std::to_string(chisqPT - chisqMin);
			// if ((chisqPT - chisqMin) > 1.99229){ // this is the crit chisq for the null pt.
			// 	// num_excluded_v[inj_idx]++;
			// }
		}//End universe search

		TFile outfile("dataTestOutFile.root", "RECREATE");
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
		// for (int inj_idx =0 ; inj_idx<4;inj_idx++){
		// 	gStyle->SetOptStat(0);
		// 	std::string name = "BestFit_"+std::to_string(inj_idx);
		// 	TCanvas can("canvas","canvas",1600,1200);
		// 	TH2D bestFits(name.c_str(),name.c_str(),25,0,25,25,0,25);
		// 	for (int expi=0; expi < nFakeExp;expi++){
		// 		bestFits.Fill(bestFitsinPt_v[inj_idx][expi],bestFitmiPt_v[inj_idx][expi]);
		// 	}
		// 	bestFits.Draw("COLZ");
		// 	std::vector<int> sidx_v = {7.5, 16.5, 19.5, 23.5};
		// 	std::vector<int> midx_v = {14.5,14.5, 14.5, 14.5};
		// 	g.SetPoint(0,sidx_v[inj_idx],midx_v[inj_idx]);
		// 	g.Draw("PSAME");
		// 	name = "BestFit_"+std::to_string(inj_idx)+".png";
		// 	can.SaveAs(name.c_str());
		// }
		//
		// std::cout <<"Ending Injection Point Test, Results:\n";
		// for (int inj_idx = 0; inj_idx < injectionPtMNu_v.size(); inj_idx++){
		// 	std::cout << "sin22: " << injectionPtSin22th_v[inj_idx] << " m2: " << injectionPtMNu_v[inj_idx] << " " << num_excluded_v[inj_idx] << " " << num_excluded_v[inj_idx]/nFakeExp << "\n";
		// }
	}//End non-gen if statement
	std::cout << "\n";
	std::cout << "End of Script \n";
	return 0;
}
