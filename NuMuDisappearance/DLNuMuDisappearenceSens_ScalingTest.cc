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

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;
std::string ZeroPadNumber(int num,int digits=3);
std::string ZeroPadNumber(int num, int digits)
{
	std::stringstream ss;

	// the number is converted to string with the help of stringstream
	ss << num;
	std::string ret;
	ss >> ret;

	// Append zero chars
	int str_length = ret.length();
	for (int i = 0; i < digits - str_length; i++)
		ret = "0" + ret;
	return ret;
}

int main(int argc, char* argv[]){
	std::string xml = "/uboone/app/users/jmills/whipping_star/xml/numudisappearance.xml";
	int iarg = 0;
	opterr=1;
	int index;
	bool gen = false;
	// const double scaleFactor(1.0);
	// const double scaleFactor(2.0);
	const double scaleFactor((190454.0/4848.0));
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
	const int nFakeExp(1000);
	const int nBins(19);
	double  mnu, sin22th, ue, um,_mnu, _sin22th, _ue, _um, chisqTest, chisqMin, chisqPT;
	int count;
	std::vector<float> fakeData;
	fakeData.resize(nBins);
	SBNspec truespec_first, appSpec_first, truespec_second, appSpec_second, testSpec;
	std::array<double,nBins>  a_pred;

	std::string s_fracDetSysMatrix = "/uboone/app/users/jmills/whipping_star/NuMuDisappearance/fracdetvar_19bins.txt";
	//PART 1: precompute all of our sin and sin2 amplitudes so we don't need to later
	if(gen){
		// use true if we need to regenerate a covariance matrix
		if(true){
		    SBNcovariance _covar(xml);
		   _covar.FormCovarianceMatrix(tag);
		   _covar.PrintMatricies(tag);
		   _covar.frac_covariance.Print();
		}
		// start with a null model
		NeutrinoModel nullModel(0, 0, 0);
		SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
		SBNspec bkg = bkgo->spec_central_value;
		//SBNspec bkg("DL.SBNspec.root",xml);
		// set full osc to 0 for bkg spectrum
		// bkg.Scale("fullosc",0.0);
		if (scaleFactor == 1.0){
			bkg.WriteOut(tag+"_Bkg");
		}
		else{
			bkg.WriteOut(tag+"_40x_Bkg");
		}


		float mnu;
		// for(int mi = 13; mi < 25; mi++){
		for(int mi = 0; mi < 13; mi++){

			mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

			//Model: mnu, ue4, um4
			//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
			//see davio's thesis eqs 4.14-16
			//we're precomputing, so we don't really care about the u's set to sin2=1
			// m41, Ue, Um
			// NeutrinoModel testModel(mnu, sqrt(.5), sqrt(.5)); //Nue App Maximized
			NeutrinoModel testModel(mnu, sqrt(.5), sqrt(.5)); //NuMu Disap Maximized

			// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
			SBNgenerate * gen = new SBNgenerate(xml,testModel);

			// Write them to files
			gen->WritePrecomputedOscSpecs(tag);
		}
		return 1;
	}

	//PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
	if(!gen){
		// Load up the necesary bits and store them in vectors on the stack **
		std::string endStr;
		if (scaleFactor == 1.0){
			endStr = "_Bkg.SBNspec.root";
		}
		else{
			endStr = "_40x_Bkg.SBNspec.root";
		}
		SBNspec cvSpec(tag+endStr,xml);
		cvSpec.CollapseVector();
		//Scale all by ScaleFactor
		cvSpec.ScaleAll(scaleFactor);
		// cvSpec.Scale("fullosc",0.0);
		// std::cout << "\n\n";
		// cvSpec.PrintCollapsedVector();
		// std::cout << "\n\n";
		// assert (1==2);
		// Setup Root file for output as well:


		std::array<SBNspec,dm2_grdpts> a_sinsqSpec;
		for(int mi = 0; mi < dm2_grdpts; mi++){
			mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			std::stringstream stream;
			stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
			std::string infile = "DL_SINSQ_dm_"+stream.str()+".SBNspec.root";
			// if (scaleFactor != 1.0){
			// 	infile = "DL_SINSQ_dm_40x_"+stream.str()+".SBNspec.root";
			// }
			auto inspec = SBNspec(infile,xml);
			inspec.ScaleAll(scaleFactor);
			// inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
			inspec.CollapseVector();
			a_sinsqSpec[mi] = inspec;
		}

		//Bring in our covariance matrix!
		// Stats only
		SBNchi uboone_chi_statsonly(cvSpec,true);

		// Stats + sys
		// Load up cov matrix and add in detector variation component	**
		TFile * fsys = new TFile("DL.SBNcovar.root","read");
		TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
		TMatrixD * covFracSys_collapsed = (TMatrixD*)fsys->Get("collapsed_frac_covariance");
		// need to add in detvar still...
		std::ifstream file;
		double ddummy;
		file.open(s_fracDetSysMatrix);
		if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << s_fracDetSysMatrix << std::endl;
		for(short i = 0; i < nBins; i++){
			for(short j = 0; j < nBins; j++){
				file >> ddummy;
				(*covFracSys_collapsed)(i,j) += ddummy;
			}
		}
		file.close();
		std::cout<<"Loaded Covariance"<<std::endl;


		TFile outfile("OutfileScalingTest_40x.root", "RECREATE");
		TH1F cv_h("cv_h","cv_h",19,250,1200);
		cv_h.SetLineColor(1);

		for (int bin=0; bin<cvSpec.collapsed_vector.size(); bin++){
			cv_h.SetBinContent(bin+1,cvSpec.collapsed_vector[bin]);
		}
		cv_h.Write();
		TDirectory *dumpDirectory = outfile.mkdir("Dump");
		// TDirectory *allOscDirectory = outfile.mkdir("Oscs");
		TDirectory *topDirectory = outfile.mkdir("GridSearch");
		// // Load up oscillation model
		std::cout << "NUMU DISAPPEARANCE" << std::endl;

		float mnu_base, um_base, ue_base, sin22th_first;
		// first get the central grid point we are testing

		// Setting up a way to call the script piecemeal.
		int miBaseStart = 0;
		int miBaseLimit = dm2_grdpts;
		// if (doOneMass){
		// 	miBaseStart = miToDo;
		// 	miBaseLimit = miToDo + 1;
		// }
		int sin22BaseStart = 0;
		int sin22BaseLimit = sin22th_grdpts;
		// if (doOneSin22){
		// 	sin22BaseStart = sin22ToDo;
		// 	sin22BaseLimit = sin22ToDo + 1;
		// }
		TDirectory *thisDirectory[2];
		for(int mi_base = miBaseStart; mi_base <miBaseLimit; mi_base++){
			for(int sin22thi_base = sin22BaseStart; sin22thi_base <sin22BaseLimit; sin22thi_base++){
				if ((((mi_base == 12) && (sin22thi_base == 20)) || ((mi_base == 19) && (sin22thi_base == 18))) == false){
					continue;
				}

				int nbinsChis = 50;
				double upperChiBound = 50;
				double upperdChiBound = 50;
				double lowerdChiBound = upperdChiBound/nbinsChis*-1.0;

				if (scaleFactor != 1.0){
					upperChiBound  = 200;
					upperdChiBound = 100;
					lowerdChiBound = upperdChiBound/nbinsChis*-1.0;
				}
				dumpDirectory->cd();
				TH1F chi2PT_dist_h("chi2PT_dist_h","chi2PT_dist_h",50,0,upperChiBound);
				TH1F chi2MIN_dist_h("chi2MIN_dist_h","chi2MIN_dist_h",50,0,upperChiBound);
				TH1F dchi2_dist_h("dchi2_dist_h","dchi2_dist_h",51,lowerdChiBound,upperdChiBound);
				std::string dirName = "massIdx_"+ZeroPadNumber(mi_base,2)+"_sinIdx_"+ZeroPadNumber(sin22thi_base,2);
				topDirectory->cd();
				int thisDirIdx;
				if ((mi_base == 12) && (sin22thi_base == 20)) {
					thisDirIdx = 0;
				}
				else{
					thisDirIdx = 1;
				}
				std::cout << thisDirIdx << " This DIR IDX " << dirName << "\n";
				thisDirectory[thisDirIdx] = topDirectory->mkdir(dirName.c_str());
				dumpDirectory->cd();

				sin22th_first = pow(10.,((sin22thi_base+0.5)/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
				mnu_base = pow(10.,((mi_base+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
				std::cout << "Running for Sin and Mass:" << sin22th_first << " " << mnu_base << "\n";
				// For NuMu Disappearance:
				um_base = 1.0 ; //THIS VAUE IS ALSO DUMMY
				ue_base = 1.0 ; // Dummy Value to satisfy code
				// current test model
				std::cout << "NU Base Model: m41^2:" << pow(mnu_base,2) << " ue:" << ue_base << " sin2:" << sin22th_first << std::endl;

				// NuMu Disappearance Math described by equations 4.20-4.22 in Davio's thesis.
				truespec_first = cvSpec;
				truespec_first.Scale("fullosc",0.0);
				std::cout<<"cv spec"<<std::endl;
				std::cout << "Before and After Spectra \n";
				truespec_first.PrintCollapsedVector();
				appSpec_first = a_sinsqSpec[mi_base];
				// appSpec_first.Scale("fullosc",sin22th_first);
				appSpec_first.Scale("bnb",-1.0*sin22th_first);
				appSpec_first.Scale("nue",0.0);
				appSpec_first.Scale("ccpi0",-1.0*sin22th_first);
				appSpec_first.Scale("ncpi0",-1.0*sin22th_first);
				appSpec_first.Scale("ext",0.0);
				truespec_first.Add(&appSpec_first);
				std::cout<<"oscillated spec"<<std::endl;
				truespec_first.PrintCollapsedVector();
				dumpDirectory->cd();
				TH1F osc_cv_h("osc_cv_h","osc_cv_h",19,250,1200);
				for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
					osc_cv_h.SetBinContent(bin+1,truespec_first.collapsed_vector[bin]);
				}
				dumpDirectory->cd();
				TH2D nearestGridPts_h("BestFit_Pts","BestFit_Pts",25,0,25,25,0,25);
				thisDirectory[thisDirIdx]->cd();
				osc_cv_h.Write();
				dumpDirectory->cd();
				if (doOneSin22){
					if ((mi_base < miToDo) || (mi_base >= miToDo+1) || (sin22thi_base < sin22ToDo) || (sin22thi_base >= sin22ToDo+1)) {
						continue;
					}
				}
				else if (doOneMass){
					if ((mi_base < miToDo) || (mi_base >= miToDo+1)) {
						continue;
					}
				}
				std::cout << "\n\n\n\n\n\n\n\nDoing This one\n\n";

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
					fakeData = TrueChi.GeneratePseudoExperiment();	// [fakedata] = vector of floats
					// fakeData = {61, 250, 302, 514, 691, 947, 846, 752, 475, 437, 278, 142, 57, 177, 108, 38, 0, 28, 43};

					// remove zero bins, they cause nans in the chi2 calc
					// for(int i = 0; i < nBins; i++){
					// 	if (fakeData[i] ==0) fakeData[i]=0.00000001;
					// }
					// print out new oscillated Spectrum
					std::cout<<"osc spec"<<std::endl;
					truespec_first.PrintCollapsedVector();
					std::cout<<"new uni spec: ";
					for(int i = 0; i < nBins; i++){
						std::cout<<fakeData[i]<<" ";
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
					double chisqMin = 99999;
					double chisqPT  = 99999;
					int count=0;
					double sin22th_second;
					int min_sin22_idx=-1;
					int min_mi_idx=-1;
					for(int mi_second = 0; mi_second <dm2_grdpts; mi_second++){
						for(int sin22thi_second = 0; sin22thi_second <sin22th_grdpts; sin22thi_second++){
							TMatrixD cov(nBins,nBins), invcov(nBins,nBins), covCopy(nBins,nBins), isDiag(nBins,nBins);
							sin22th_second = pow(10.,((sin22thi_second+0.5)/float(sin22th_grdpts)*TMath::Log10(sin22th_hibound/sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
							// for(int i = 0; i < nBins; i++){
							// 	a_pred[i] = cvSpec.collapsed_vector[i] - _sin22th*a_sinsqSpec[mi].collapsed_vector[i];
							// }

							//
							truespec_second = cvSpec;
							truespec_second.Scale("fullosc",0.0);
							appSpec_second = a_sinsqSpec[mi_second];
							// appSpec_second.Scale("fullosc",sin22th_second);
							appSpec_second.Scale("bnb",-1.0*sin22th_second);
							appSpec_second.Scale("nue",0.0);
							appSpec_second.Scale("ccpi0",-1.0*sin22th_second);
							appSpec_second.Scale("ncpi0",-1.0*sin22th_second);
							appSpec_second.Scale("ext",0.0);
							truespec_second.Add(&appSpec_second);



							// haven't touched shape only yet, but leaving in from davio for use later.
							if(shapeonly){

								double covIntegral(0.0), predIntegral(0.0), obsIntegral(0.0),fnorm;
								for(int i = 0; i < nBins; i++){
									predIntegral+=truespec_second.collapsed_vector[i];
									obsIntegral+=fakeData[i];
									for(int j = 0; j < nBins; j++){
										covIntegral += (*covFracSys_collapsed)[i][j]*truespec_second.collapsed_vector[i]*truespec_second.collapsed_vector[j];
									}
								}
								fnorm = covIntegral/pow(predIntegral,2);
								for(int i = 0; i < nBins; i++){
									truespec_second.collapsed_vector[i] *= (obsIntegral/predIntegral); // normalize prediction
								}
								for(int i = 0; i < nBins; i++){
									for(int j = 0; j < nBins; j++){
										cov[i][j] = ((*covFracSys_collapsed)[i][j]-fnorm)*truespec_second.collapsed_vector[i]*truespec_second.collapsed_vector[j];
										if(i==j){
											if (truespec_second.collapsed_vector[i] >0 && fakeData[i] >0){
												cov[i][j] += 3.0 / (1.0/fakeData[i] + 2.0/truespec_second.collapsed_vector[i]);
												// poisson version
												// cov = (M-mu)**2 / (2*(mu - M + M*np.log(M/mu)))
												// cov[i][j] += pow(fakeData[i]-truespec_second.collapsed_vector[i],2) / (2.0*(truespec_second.collapsed_vector[i]- fakeData[i]+ fakeData[i]*log(fakeData[i]/truespec_second.collapsed_vector[i])));
											}
											else {
												cov[i][j] += truespec_second.collapsed_vector[i]/2.0;
											}
											// Old way
											// cov[i][j] += truespec_second[i];
										}
									}
								}
							}//End of shape only
							// Shape+Rate
							else{
								for(int i = 0; i < nBins; i++){
									for(int j = 0; j < nBins; j++){
											// remove some nans
											if(truespec_second.collapsed_vector[i] > 0 && truespec_second.collapsed_vector[j] >0 ){
												cov[i][j] += (*covFracSys_collapsed)[i][j]*truespec_second.collapsed_vector[i]*truespec_second.collapsed_vector[j];
											}

											if(i==j){
												// cov[i][j] += truespec_second.collapsed_vector[i]; // Davio includes this
												// add in stat errors start with CNP for "data" errors
												if (truespec_second.collapsed_vector[i] >0 && fakeData[i] >0){
													cov[i][j] += 3.0 / (1.0/fakeData[i] + 2.0/truespec_second.collapsed_vector[i]);
													// poisson version
													// cov = (M-mu)**2 / (2*(mu - M + M*np.log(M/mu)))
													// cov[i][j] += pow(fakeData[i]-truespec_second.collapsed_vector[i],2) / (2.0*(truespec_second.collapsed_vector[i]- fakeData[i]+ fakeData[i]*log(fakeData[i]/truespec_second.collapsed_vector[i])));
												}
												else {
													cov[i][j] += truespec_second.collapsed_vector[i]/2.0;
												}
												// mc stat error
												// cov[i][j] +=stkerr[i]*stkerr[i];
											}
										}
								}
							} //end of shape+rate
							// inv cov for chi2calc
							covCopy = cov;
							invcov = covCopy.Invert();
							isDiag.Mult(invcov,cov);

							chisqTest = 0;
							for(int i = 0; i < nBins; i++){
								for(int j = 0; j < nBins; j++){
									// (obsi-predi)*(invcov)*(obsj-predj)
									double addVal = (fakeData[i] - truespec_second.collapsed_vector[i])*invcov[i][j]*(fakeData[j] - truespec_second.collapsed_vector[j]);
									chisqTest += addVal;
								}
							}
							if ((mi_second == mi_base) && (sin22thi_base == sin22thi_second)){
								chisqPT = chisqTest;
							}
							if (chisqTest < chisqMin){
								chisqMin = chisqTest;
								min_sin22_idx = sin22thi_second;
								min_mi_idx = mi_second;
							}
						}//End second sin22 search
					}//End second m41 search
					std::cout << chisqPT << " " << chisqMin << " " << chisqPT - chisqMin << " Chisquares\n";
					// Fill Hists:
					if (chisqPT >= upperChiBound){
						chi2PT_dist_h.Fill(upperChiBound-.1);
					}
					else{
						chi2PT_dist_h.Fill(chisqPT);
					}
					if (chisqMin >= upperChiBound){
						chi2MIN_dist_h.Fill(upperChiBound-.1);
					}
					else{
						chi2MIN_dist_h.Fill(chisqMin);
					}
					if ((chisqPT - chisqMin) >= upperdChiBound){
						dchi2_dist_h.Fill(upperdChiBound-.1);
					}
					else if ((chisqPT - chisqMin) < lowerdChiBound){
						dchi2_dist_h.Fill(lowerdChiBound);
					}
					else{
						dchi2_dist_h.Fill(chisqPT - chisqMin);
					}
					std::cout << "Min Sin22: "<< min_sin22_idx << " Min Mi:" << min_mi_idx << "\n";
					nearestGridPts_h.Fill(min_sin22_idx,min_mi_idx);
					std::string histName = "uni_"+ZeroPadNumber(expi,4)+"_dChi2_"+std::to_string(chisqPT - chisqMin);
					TH1F uni_h(histName.c_str(),histName.c_str(),19,250,1200);
					uni_h.SetLineColor(2);
					for (int bin=0; bin<fakeData.size(); bin++){
						uni_h.SetBinContent(bin+1,fakeData[bin]);
					}
					thisDirectory[thisDirIdx]->cd();
					uni_h.Write();
					dumpDirectory->cd();
				}//End universe search
				thisDirectory[thisDirIdx]->cd();
				nearestGridPts_h.Write();
				chi2PT_dist_h.Write();
				chi2MIN_dist_h.Write();
				dchi2_dist_h.Write();
				thisDirectory[thisDirIdx]->Write();
				// outfile.cd();
				outfile.Write();

			}
		}
	}
	std::cout << "End of Script \n";
	return 0;
}
