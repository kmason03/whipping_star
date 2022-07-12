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
#include "TPaveText.h"

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
	const double scaleFactor(1.0);
	// const double scaleFactor(2.0);
	// const double scaleFactor((190454.0/4848.0));
	bool numudis = false;
	bool nueapp = true;
	bool combined = false;
	int mass_start = -1;
	bool shapeonly = true;
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
	{"statsOnly", no_argument,0,'st'},
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
			case 'st':
				statsOnly = true;
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
	const int nFakeExp(1000);
	const int nBins(19);
	double  mnu, sin22th, ue, um,_mnu, _sin22th, _ue, _um;
	int count;
	std::vector<double> fakeData;
	fakeData.resize(nBins);
	SBNspec truespec_first, appSpec_first, truespec_second, appSpec_second, testSpec;
	std::array<double,nBins>  a_pred;

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	if(!gen){
		// save output to text files
		std::ofstream chifile;
		std::string chiFileName = "rterms.txt";
		if (statsOnly){
			if (doOneMass){
				if (!shapeonly){
					chiFileName = "rterms_statsonly_1x_pt_"+ZeroPadNumber(miToDo)+".txt";
				}
				else{
					chiFileName = "rterms_statsonly_shapeonly_1x_pt_"+ZeroPadNumber(miToDo)+".txt";
				}
			}
			if (doOneSin22){
				if (!shapeonly){
					chiFileName = "rterms_statsonly_pt_"+ZeroPadNumber(miToDo)+"_"+ZeroPadNumber(sin22ToDo)+".txt";
				}
				else{
					chiFileName = "rterms_statsonly_shapeonly_pt_"+ZeroPadNumber(miToDo)+"_"+ZeroPadNumber(sin22ToDo)+".txt";
				}
			}
		}
		else{
			if (doOneMass){
				if (!shapeonly){
					chiFileName = "rterms_1x_pt_"+ZeroPadNumber(miToDo)+".txt";
				}
				else{
					chiFileName = "rterms_shapeonly_1x_pt_"+ZeroPadNumber(miToDo)+".txt";
				}
			}
			if (doOneSin22){
				if (!shapeonly){
					chiFileName = "rterms_pt_"+ZeroPadNumber(miToDo)+"_"+ZeroPadNumber(sin22ToDo)+".txt";
				}
				else{
					chiFileName = "rterms_shapeonly_pt_"+ZeroPadNumber(miToDo)+"_"+ZeroPadNumber(sin22ToDo)+".txt";
				}
			}
		}

		dumpCoordFile(dm2_grdpts, dm2_lowbound, dm2_hibound, sin22th_grdpts, sin22th_lowbound);

		chifile.open(chiFileName.c_str(), std::ios_base::app);
		// Load up the necesary bits and store them in vectors on the stack **
		std::string endStr = "_Bkg.SBNspec.root";
		// if (scaleFactor == 1.0){
		// 	endStr = "_Bkg.SBNspec.root";
		// }
		// else{
		// 	endStr = "_1x_Bkg.SBNspec.root";
		// }
		SBNspec cvSpec(tag+endStr,xml);
		cvSpec.RemoveMCError();
		cvSpec.CollapseVector();
		//Scale all by ScaleFactor
		cvSpec.ScaleAll(scaleFactor);



		std::vector<SBNspec> a_sinsqSpec;
		getPrecompSpec(a_sinsqSpec, dm2_grdpts, dm2_hibound, dm2_lowbound, xml, scaleFactor);
		// assert (1==2);


		//Bring in our covariance matrix!
		// Stats + sys
		// Load up cov matrix and add in detector variation component	**
		std::string covInFile;
		if (statsOnly) {covInFile = "JOSH_1m1p_StatOnly.SBNcovar.root";}
		else {covInFile = "JOSH_1m1p.SBNcovar.root";}
		TFile * fsys = new TFile(covInFile.c_str(),"read");
		TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
		TMatrixD * covFracSys_collapsed = (TMatrixD*)fsys->Get("frac_covariance_collapsed");
		// need to add in detvar still...
		std::cout<<"Loaded Covariance"<<std::endl;
		if (false){
			gStyle->SetOptStat(0);
			TH2F covar_h("covar_h","covar_h",19,250,1200,19,250,1200);
			for (int a=0; a<19;a++){
				for (int b=0; b<19;b++){
					covar_h.SetBinContent(a+1,b+1,(*covFracSys_collapsed)[a][b]);
				}
			}
			//covar_h.SetMaximum(0.07);
			TCanvas tmpcan("tmpcan","tmpcan",1200,1000);
			tmpcan.SetRightMargin(0.12);
			tmpcan.SetLeftMargin(0.12);
			covar_h.SetTitle("Stats-Only  Covariance Matrix");
			covar_h.SetXTitle("Neutrino Energy");
			covar_h.SetYTitle("Neutrino Energy");
			covar_h.Draw("COLZ");
			tmpcan.SaveAs("FractionalCovarMat_Total.png");
			return 0;
		}



		std::string outfile_name;
		if (statsOnly) {outfile_name = "Outfile_statsonly_"+ZeroPadNumber(miToDo)+".root";}
		else {outfile_name = "Outfile_statsonly_"+ZeroPadNumber(miToDo)+".root";}

		TFile outfile(outfile_name.c_str(), "RECREATE");
		TH1F cv_h("cv_h","cv_h",19,250,1200);
		cv_h.SetLineColor(1);

		for (int bin=0; bin<cvSpec.collapsed_vector.size(); bin++){
			cv_h.SetBinContent(bin+1,cvSpec.collapsed_vector[bin]);
		}
		cv_h.Write();

		TH1F chi_dist_h("chi_dist_h","chi_dist_h",25,0,200);

		TDirectory *dumpDirectory = outfile.mkdir("Dump");
		TDirectory *topDirectory = outfile.mkdir("GridSearch");

		// // Load up oscillation model
		std::cout << "NUMU DISAPPEARANCE" << std::endl;

		float mnu_base, um_base, ue_base, sin22th_first;
		// first get the central grid point we are testing

		// Setting up a way to call the script piecemeal.
		int miBaseStart = 0;
		int miBaseLimit = dm2_grdpts;
		int sin22BaseStart = 0;
		int sin22BaseLimit = sin22th_grdpts;

		TDirectory *thisDirectory[(miBaseLimit-miBaseStart)*(sin22BaseLimit-sin22BaseStart)];
		std::cout << miBaseStart << " " << miBaseLimit << " "  << sin22BaseStart << " " << sin22BaseLimit << "\n";


		std::vector<std::vector<TMatrixD> > fracCovMats_grid = genCovMats(cvSpec,
																								a_sinsqSpec,
																								*covFracSys,
																								dm2_grdpts,
																								sin22th_grdpts,
																								sin22th_lowbound,
																								nBins
																							);


		for(int mi_base = miBaseStart; mi_base <miBaseLimit; mi_base++){
			for(int sin22thi_base = sin22BaseStart; sin22thi_base <sin22BaseLimit; sin22thi_base++){
				// if ((sin22thi_base < 12) || (sin22thi_base >18)){
				// 	continue;
				// }
				std::string dirName = "massIdx_"+ZeroPadNumber(mi_base,2)+"_sinIdx_"+ZeroPadNumber(sin22thi_base,2);
				topDirectory->cd();
				int thisDirIdx = (mi_base-miBaseStart)*(miBaseLimit-miBaseStart) + (sin22thi_base-sin22BaseStart);
				std::cout << thisDirIdx << " Dir Idx\n";
				thisDirectory[thisDirIdx] = topDirectory->mkdir(dirName.c_str());
				std::cout << " Dir Idx AGAIN\n";
				dumpDirectory->cd();


					//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
				sin22th_first = pow(10.,((sin22thi_base+0.5)/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
				mnu_base = pow(10.,((mi_base+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
				std::cout << "Running for Sin and Mass:" << sin22th_first << " " << mnu_base << "\n";
				// For Nue Appearance:
				// um_base = .135;	// global best fit
				// ue_base = pow(sin22th_first/float(4*um_base*um_base),.5);
				// For NuMu Disappearance:
				um_base = 1.0 ; //THIS VAUE IS ALSO DUMMY
				ue_base = 1.0 ; // Dummy Value to satisfy code
				// current test model
				std::cout << "NU Base Model: m41^2:" << pow(mnu_base,2) << " ue:" << ue_base << " sin2:" << sin22th_first << std::endl;

				// NuMu Disappearance Math described by equations 4.20-4.22 in Davio's thesis.
				truespec_first = cvSpec;
				std::cout<<"cv spec"<<std::endl;
				cvSpec.PrintCollapsedVector();
				std::cout << mi_base << " " << a_sinsqSpec.size() << "\n";
				scaleDisappearedSpec(truespec_first, sbn::SBNspec(a_sinsqSpec[mi_base]), sin22th_first);
				std::cout<<"disappeared spec"<<std::endl;
				truespec_first.PrintCollapsedVector();
				TH1F osc_cv_h("osc_cv_h","osc_cv_h",19,250,1200);
				for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
					osc_cv_h.SetBinContent(bin+1,truespec_first.collapsed_vector[bin]);
				}
				thisDirectory[thisDirIdx]->cd();
				TH2D nearestGridPts_h("BestFit_Pts","BestFit_Pts",sin22th_grdpts,0,sin22th_grdpts,dm2_grdpts,0,dm2_grdpts);
				osc_cv_h.Write();
				dumpDirectory->cd();
				if (doOneSin22){
					if ((mi_base < miToDo) || (mi_base >= miToDo+1) || (sin22thi_base < sin22ToDo) || (sin22thi_base >= sin22ToDo+1)) {
						continue;
					}
				}
				else{
					if ((mi_base < miToDo) || (mi_base >= miToDo+1)) {
						continue;
					}
				}
				std::cout << "\n\n\n\n\nDoing This one\n\n";

				//  I'll be generating a new universe around the oscillated spectrum
				// std::cout<<"truespec size: "<<truespec_first.collapsed_vector.size()<<std::endl;
				// std::cout<<"matrix size:   "<<covFracSys->GetNrows()<<std::endl;

				 // I'll be generating a new universe around the oscillated spectrum

				// SBNchi TrueChi(truespec_first,*covFracSys);
				// TrueChi.ReloadCoreSpectrum(&truespec_first);
				// TMatrixD covFracSys_base = getZeroedMat(nBins);
				// TrueChi.FillCollapsedFractionalMatrix(&covFracSys_base);
				SBNchi TrueChi(truespec_first,*covFracSys);
				TMatrixD covFracSys_test = getZeroedMat(nBins);
				TrueChi.FillCollapsedFractionalMatrix(&covFracSys_test);
				if (false){
					gStyle->SetOptStat(0);
					TH2F covar_h("covar_h","covar_h",19,250,1200,19,250,1200);
					for (int a=0; a<19;a++){
						for (int b=0; b<19;b++){
							covar_h.SetBinContent(a+1,b+1,(fracCovMats_grid[mi_base][sin22thi_base])[a][b]);
							covar_h.SetBinContent(a+1,b+1,covFracSys_test[a][b]);

						}
					}
					//covar_h.SetMaximum(0.07);
					TCanvas tmpcan("tmpcan","tmpcan",1200,1000);
					covar_h.Draw("COLZ");
					tmpcan.SaveAs("FractionalCovarMat_Thrown.png");
					osc_cv_h.Draw();
					tmpcan.SaveAs("F_Osc.png");
					cv_h.Draw();
					tmpcan.SaveAs("F_CV.png");
				}
				TrueChi.pseudo_from_collapsed = true;
				TrueChi.GeneratePseudoExperiment();		// get the motor running with initial cholosky decomposition

				TH1F fakedata_h("fakedata_h","fakedata_h",19,250,1200);

				// Across several fake experiments:
				for(int expi = 0; expi < nFakeExp; expi++){
					std::cout << expi << " Experiment\n";
					std::vector<float> fakeData_f = TrueChi.GeneratePseudoExperiment();	// [fakedata] = vector of floats
					for (int bini = 0;bini<fakeData.size();bini++){
						fakeData[bini] = (double)fakeData_f[bini];
						fakedata_h.SetBinContent(bini+1,fakeData_f[bini]);
					}

					// print out new oscillated Spectrum
					std::cout<<"osc spec"<<std::endl;
					truespec_first.PrintCollapsedVector();
					std::cout<<"new uni spec: ";
					for(int i = 0; i < nBins; i++){
						std::cout<<fakeData[i]<<" ";
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
					TH1F bestfitspec_h("bestfitspec_h","bestfitspec_h",19,250,1200);
					for(int mi_second = 0; mi_second <dm2_grdpts; mi_second++){
						for(int sin22thi_second = 0; sin22thi_second <sin22th_grdpts; sin22thi_second++){
							TMatrixD cov(nBins,nBins), invcov(nBins,nBins), covCopy(nBins,nBins), isDiag(nBins,nBins);
							sin22th_second = pow(10.,((sin22thi_second+0.5)/float(sin22th_grdpts)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
							mnu_second = pow(10.,((mi_second+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

							truespec_second = cvSpec;
							if ((false) && (sin22th_second > 0.9)){
								std::cout <<"---------------------------------------------------\n";
								std::cout << mi_second << " " << sin22th_second << "\n";
							}
							scaleDisappearedSpec(truespec_second, sbn::SBNspec(a_sinsqSpec[mi_second]), sin22th_second, false);
							// SBNchi TrueChi_Second(truespec_second,*covFracSys);
							// // TrueChi_Second.ReloadCoreSpectrum(&truespec_second);
							// TMatrixD covFracSys_second = getZeroedMat(nBins);
							// TrueChi_Second.FillCollapsedFractionalMatrix(&covFracSys_second);

							// haven't touched shape only yet, but leaving in from davio for use later.
							if(shapeonly){
								scaleCovShapeOnly(cov, fracCovMats_grid[mi_second][sin22thi_second], fakeData, truespec_second.collapsed_vector, nBins); // cov, covfraccollapsed, obs, exp(gets modified), nbins
							}//End of shape only
							// Shape+Rate
							else{
								scaleCovShapeRate(cov, fracCovMats_grid[mi_second][sin22thi_second], fakeData, truespec_second.collapsed_vector, nBins); // cov, covfraccollapsed, obs, exp, nbins
							} //end of shape+rate

							if ((false) && (expi == 0) && (mi_second == 12) && (sin22thi_second == 0)){
								gStyle->SetOptStat(0);
								TH2F covar_h("covar_h","covar_h",19,250,1200,19,250,1200);
								for (int a=0; a<19;a++){
									for (int b=0; b<19;b++){
										covar_h.SetBinContent(a+1,b+1,(fracCovMats_grid[mi_second][sin22thi_second])[a][b]);
									}
								}
								//covar_h.SetMaximum(0.07);
								TCanvas tmpcan("tmpcan","tmpcan",1200,1000);
								covar_h.Draw("COLZ");
								tmpcan.SaveAs("FractionalCovarMat_Used.png");
								// assert(1==2);
							}

							covCopy = cov;
							invcov = covCopy.Invert();
							isDiag.Mult(invcov,cov);
							calcRTerm(detTermTest, chisqTest, fakeData, truespec_second.collapsed_vector, invcov, cov, nBins); // obs, exp, inv, nbins
							halfRTermTest = detTermTest + chisqTest;
							if ((mi_second == mi_base) && (sin22thi_base == sin22thi_second)){
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
								for (int bin=0; bin<bestfitspec_h.GetNbinsX(); bin++){
									bestfitspec_h.SetBinContent(bin+1,truespec_second.collapsed_vector[bin]);
								}
							}
						}//End second sin22 search
					}//End second m41 search

					if (false){
						TCanvas tmpcan("tmpcan","tmpcan",1200,1000);
						cv_h.SetLineColor(kBlack);
						cv_h.SetMaximum(450);
						cv_h.Draw("");
						osc_cv_h.SetLineColor(kBlue);
						osc_cv_h.Draw("SAME");
						fakedata_h.SetLineColor(kRed);
						fakedata_h.Draw("SAME");

						TPaveText paveText(400,30,900,155,"NB");
				    // col = ROOT.TColor()
				    // col.SetRGB(0,0,0)
				    // paveText.SetFillColorAlpha(col.GetNumber(),1.0)
				    paveText.SetTextAlign(13);
						paveText.AddText(("PT Mnu^2: " + std::to_string(pow(mnu_base,2)) + " Sin^2: "+std::to_string(sin22th_first)).c_str());
						paveText.AddText(("BF Params: " + std::to_string(pow(best_mnu_second,2)) + " Sin^2: "+std::to_string(best_sin22th_second)).c_str());
				    paveText.AddText(("PT Det: " + std::to_string(detTermPT)).c_str());
						paveText.AddText(("BF Det: " + std::to_string(detTermMin)).c_str());
						paveText.AddText(("PT Chi2:" + std::to_string(chisqPT)).c_str());
						paveText.AddText(("BF Chi2:" + std::to_string(chisqMin)).c_str());
						paveText.AddText(("PT Tot: " + std::to_string(detTermPT+chisqPT)).c_str());
						paveText.AddText(("BF Tot: " + std::to_string(detTermMin+chisqMin)).c_str());




						TLegend legend(0.7,0.70,0.90,0.90);
						legend.AddEntry(&cv_h,					"Null Osc   ");
						legend.AddEntry(&osc_cv_h,			"PT Osc     ");
						legend.AddEntry(&fakedata_h,		"Thrown Uni ");

						legend.Draw("SAME");
						tmpcan.SaveAs(("unimatches/FakeUni_" +ZeroPadNumber(sin22thi_base)+"_"+ZeroPadNumber(expi)+"_Alone.png").c_str());
						cv_h.Draw("");
						osc_cv_h.Draw("SAME");
						fakedata_h.Draw("SAME");
						bestfitspec_h.SetLineColor(kMagenta);
						bestfitspec_h.Draw("SAME");
						legend.AddEntry(&bestfitspec_h,	"Best Osc   ");
						legend.Draw("SAME");
						paveText.Draw("SAME");
						tmpcan.SaveAs(("unimatches/FakeUni_" +ZeroPadNumber(sin22thi_base)+"_"+ZeroPadNumber(expi)+"_Best.png").c_str());

					}


					chifile << halfRTermPT - halfRTermMin << " ";
					std::cout << halfRTermPT << " " << halfRTermMin << " " << halfRTermPT - halfRTermMin << " R Terms\n";
					std::cout << detTermPT << ", " << chisqPT << ", " << detTermMin << ", " << chisqMin  << ", Term Breakdowns,";

					std::cout << "Min Sin22, "<< min_sin22_idx << ", Min Mi," << min_mi_idx << "\n";
					nearestGridPts_h.Fill(min_sin22_idx,min_mi_idx);
					std::string histName = "uni_"+ZeroPadNumber(expi,4)+"_dChi2_"+std::to_string(halfRTermPT - halfRTermMin);
					TH1F uni_h(histName.c_str(),histName.c_str(),19,250,1200);
					uni_h.SetLineColor(2);
					for (int bin=0; bin<fakeData.size(); bin++){
						uni_h.SetBinContent(bin+1,fakeData[bin]);
					}
					thisDirectory[thisDirIdx]->cd();
					// uni_h.Write();
					dumpDirectory->cd();

				}//End universe search


				if (true){
					gStyle->SetOptStat(0);
					TCanvas tmpcan("tmpcan","tmpcan",1200,1000);
					TGraph tgraph;
					nearestGridPts_h.Draw("COLZ");
					tgraph.SetPoint(0,sin22thi_base+0.5,mi_base+0.5);
					tgraph.SetMarkerSize(2);
					tgraph.SetMarkerStyle(kFullCircle);
					tgraph.SetMarkerColor(kRed);
					tgraph.Draw("SAMEP");
					tmpcan.SaveAs(("unimatches/BestFitUni_"+ZeroPadNumber(mi_base) + "_"+ZeroPadNumber(sin22thi_base)+".png").c_str());
				}


			nearestGridPts_h.Write();
			chi_dist_h.Write();
			thisDirectory[thisDirIdx]->Write();
			// outfile.cd();
			outfile.Write();

			chifile<< std::endl;


			}
		}
	}
	std::cout << "End of Script \n";
	return 0;
}
