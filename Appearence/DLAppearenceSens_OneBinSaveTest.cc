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
	std::string xml = "/uboone/app/users/kmason/whipping_star/xml/nueappearence.xml";
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
	const int nFakeExp(1000);
	const int nBins(22);
	double  mnu, sin22th, ue, um,_mnu, _sin22th, _ue, _um, chisqTest, chisqMin, chisqPT;
	int count;
	std::vector<float> fakeData;
	fakeData.resize(nBins);
	SBNspec truespec_first, appSpec_first, truespec_second, appSpec_second, testSpec;
	std::array<double,nBins>  a_pred;

	std::string s_fracDetSysMatrix = "/uboone/app/users/kmason/whipping_star/build/Appearence/fracdetvar.txt";
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
		bkg.WriteOut(tag+"_Bkg");

		float mnu;
		for(int mi = 13; mi < 25; mi++){
		// for(int mi = 0; mi < 13; mi++){

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
		// save output to text files
		std::ofstream coordfile;
		std::ofstream chifile;
		std::string chiFileName = "chis.txt";
		if (doOneMass){
			chiFileName = "chis_pt_"+ZeroPadNumber(miToDo)+".txt";
		}
		if (doOneSin22){
			chiFileName = "chis_pt_"+ZeroPadNumber(miToDo)+ZeroPadNumber(sin22ToDo)+".txt";
		}
		coordfile.open("bins.txt", std::ios_base::app);//std::ios_base::app
		chifile.open(chiFileName.c_str(), std::ios_base::app);
		coordfile.open("bins.txt", std::ios_base::app);//std::ios_base::app
		chifile.open(chiFileName.c_str(), std::ios_base::app);
		// Print binedges for easier plotting **
		for(int mi = 0; mi <= dm2_grdpts; mi++){
			mnu = pow(10.,((mi)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			coordfile << pow(mnu,2) << " ";
		}
		coordfile << std::endl;



		for(int si = 0; si <= sin22th_grdpts; si++){
			sin22th = pow(10.,(si/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
			coordfile << sin22th << " ";
		}
		coordfile << std::endl;

		// Load up the necesary bits and store them in vectors on the stack **
		SBNspec cvSpec("MassSpectraData/AppOnly/"+tag+"_Bkg.SBNspec.root",xml);
		cvSpec.Scale("fullosc",0.0);

		std::array<SBNspec,dm2_grdpts> a_sinsqSpec;
		for(int mi = 0; mi < dm2_grdpts; mi++){
			mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
			std::stringstream stream;
			stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
			std::string infile = "MassSpectraData/AppOnly/"+tag+"_SINSQ_dm_"+stream.str()+".SBNspec.root";
			auto inspec = SBNspec(infile,xml);
			// inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
			inspec.CollapseVector();
			a_sinsqSpec[mi] = inspec;
		}

		//Bring in our covariance matrix!
		// Stats only
		SBNchi uboone_chi_statsonly(cvSpec,true);

		// Stats + sys
		// Load up cov matrix and add in detector variation component	**
		TFile * fsys = new TFile("MassSpectraData/AppOnly/DL.SBNcovar.root","read");
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


		TFile outfile("Outfile.root", "RECREATE");
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
		// lets just look at nue app for now

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
		TDirectory *thisDirectory[(miBaseLimit-miBaseStart)*(sin22BaseLimit-sin22BaseStart)];
		std::cout << miBaseStart << " " << miBaseLimit << " "  << sin22BaseStart << " " << sin22BaseLimit << "\n";

		for(int mi_base = miBaseStart; mi_base <miBaseLimit; mi_base++){
			for(int sin22thi_base = sin22BaseStart; sin22thi_base <sin22BaseLimit; sin22thi_base++){

		// for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
		// 	for(int sin22thi_base = 0; sin22thi_base <sin22th_grdpts; sin22thi_base++){

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
				// Davio
				// truespec_first = cvSpec;
				// disSpec = a_sinsqSpec[pTdm2i];
				// disSpec.Scale("bnb",-sin22th);
				// truespec_first.Add(&disSpec);
				// truespec_first.CollapseVector();

				truespec_first.Scale("fullosc",0.0);
				std::cout<<"cv spec"<<std::endl;
				std::cout << "Before and After Spectra \n";
				truespec_first.PrintCollapsedVector();
				appSpec_first = a_sinsqSpec[mi_base];
				// appSpec_first.Scale("fullosc",sin22th_first);
				appSpec_first.Scale("fullosc",sin22thi_base);
				appSpec_first.Scale("bnb",0.0);
				appSpec_first.Scale("nue",0.0);
				appSpec_first.Scale("ccpi0",0.0);
				appSpec_first.Scale("ncpi0",0.0);
				appSpec_first.Scale("ext",0.0);
				truespec_first.Add(&appSpec_first);
				std::cout<<"oscillated spec"<<std::endl;
				truespec_first.PrintCollapsedVector();
				TH1F osc_cv_h("osc_cv_h","osc_cv_h",19,250,1200);
				for (int bin=0; bin<truespec_first.collapsed_vector.size(); bin++){
					osc_cv_h.SetBinContent(bin+1,truespec_first.collapsed_vector[bin]);
				}
				thisDirectory[thisDirIdx]->cd();
				TH2D nearestGridPts_h("BestFit_Pts","BestFit_Pts",25,0,25,25,0,25);
				osc_cv_h.Write();
				dumpDirectory->cd();
				if ((mi_base < miToDo) || (mi_base >= miToDo+1) || (sin22thi_base < sin22ToDo) || (sin22thi_base >= sin22ToDo+1)) {
					continue;
				}
				std::cout << "\n\n\n\n\n\n\n\nDoing This one\n\n";
				//  I'll be generating a new universe around the oscillated spectrum
				SBNchi TrueChi(truespec_first,*covFracSys);
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
							appSpec_second.Scale("fullosc",sin22thi_second);
							appSpec_second.Scale("bnb",0.0);
							appSpec_second.Scale("nue",0.0);
							appSpec_second.Scale("ccpi0",0.0);
							appSpec_second.Scale("ncpi0",0.0);
							appSpec_second.Scale("ext",0.0);
							truespec_second.Add(&appSpec_second);



							// haven't touched shape only yet, but leaving in from davio for use later.
							if(shapeonly){
								// double covIntegral(0.0), predIntegral(0.0), obsIntegral(0.0),fnorm;
								// for(int i = 0; i < nBins; i++){
								// 	obsIntegral+=a_pred[i];
								// 	predIntegral+=fakeData[i];
								// 	for(int j = 0; j < nBins; j++){
								// 		covIntegral += (*covFracSys)[i][j]*cvSpec.collapsed_vector[i]*cvSpec.collapsed_vector[j];
								// 	}
								// }
								// fnorm = covIntegral/pow(predIntegral,2);
								// for(int i = 0; i < nBins; i++){
								// 	a_pred[i] *= (obsIntegral/predIntegral); // normalize prediction
								// }
								// for(int i = 0; i < nBins; i++){
								// 	for(int j = 0; j < nBins; j++){
								// 		cov[i][j] = ((*covFracSys)[i][j]-fnorm)*a_pred[i]*a_pred[j];
								// 		if(i==j)
								// 			cov[i][j] += a_pred[i];
								// 	}
								// }
							}
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
							} //end of shape only
							// inv cov for chi2calc


							// for (int i=0;i<nBins;i++){
							// 	cov[i][i] = 5;
							// }
							covCopy = cov;
							invcov = covCopy.Invert();
							isDiag.Mult(invcov,cov);

							// std::cout << "\nCov:\n";
							// for(int i = 0; i < nBins; i++){
							// 	for(int j = 0; j < nBins; j++){
							// 		std::cout << cov[i][j] << " ";
							// 	}
							// 	std::cout << "\n";
							// }
							// std::cout << "\nInv:\n";
							// for(int i = 0; i < nBins; i++){
							// 	for(int j = 0; j < nBins; j++){
							// 		std::cout << invcov[i][j] << " ";
							// 	}
							// 	std::cout << "\n";
							// }
							// if ((mi_second == 0) && (sin22thi_second == 0)){
							// 	std::cout << "\nisDiag:\n";
							// 	for(int i = 0; i < nBins; i++){
							// 		for(int j = 0; j < nBins; j++){
							// 			double val =0;
							// 			for (int k=0;k<nBins;k++){
							// 				val = val + invcov[k][j]*cov[i][k];
							// 			}
							// 			// std::cout << val << " ";
							// 			std::cout << isDiag[i][j] << " ";
							// 		}
							// 		std::cout << "\n";
							// 	}
							// }



							chisqTest = 0;
							for(int i = 0; i < nBins; i++){
								for(int j = 0; j < nBins; j++){
									// (obsi-predi)*(invcov)*(obsj-predj)
									double addVal = (fakeData[i] - truespec_second.collapsed_vector[i])*invcov[i][j]*(fakeData[j] - truespec_second.collapsed_vector[j]);
									chisqTest += addVal;
									// Debugging why chi2 is so good for 0/0 case
									// if (((sin22thi_second == 24) && (mi_second ==15)) || ((sin22thi_second == 0) && (mi_second == 0))){
									// 	std::cout << i << " " << j << " " << fakeData[i] - truespec_second.collapsed_vector[i] << " " << invcov[i][j] << " " << fakeData[j] - truespec_second.collapsed_vector[j] << " " << addVal << " " << chisqTest<< "\n";
									// }
								}
							}
							if ((mi_second == mi_base) && (sin22thi_base == sin22thi_second)){
								chisqPT = chisqTest;
								double val = chisqPT;
								if (chisqPT > 200){
									val = 199;
								}
								chi_dist_h.Fill(val);
							}
							if (chisqTest < chisqMin){
								chisqMin = chisqTest;
								min_sin22_idx = sin22thi_second;
								min_mi_idx = mi_second;
							}
						}//End second sin22 search
					}//End second m41 search
					chifile << chisqPT - chisqMin << " ";
					std::cout << chisqPT << " " << chisqMin << " " << chisqPT - chisqMin << " Chisquares\n";
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
