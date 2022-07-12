#ifndef MILLSFUNCTIONS_H_
#define MILLSFUNCTIONS_H_

#include <string>
#include <utility>
#include "assert.h"
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
#include "TPad.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "MillsFunctions.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"
#include "SBNcovariance.h"
#include "prob.h"
#include <assert.h>

// namespace sbn{


	/**********************************************
	 *	Centralized Function Hub for J.Mills Thesis Work
	 * *******************************************/

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 3);
// Test Function for Build:
void testFunction();
// Pad String with extra zeros
std::string ZeroPadNumber(int num, int digits=3);
// Generate SBNspec and save based on input massid (int from assumed list)
void generate_spectra(int massid, std::string xml, std::string tag, double dm2_lowbound, double dm2_hibound, int dm2_grdpts);
// Overload Generate Spectra for Sterile Decay:
void generate_spectra(int massid, int Um4sq, std::string xml, std::string tag, double dm2_lowbound, double dm2_hibound, int dm2_grdpts, double Um4sq_lowbound, double Um4sq_hibound, int Um4sq_grdpts);
// Calculate the chi2 difference
void calcRTerm(double& detTerm, double& chiSQ, const std::vector<double>& obs, const std::vector<double>& exp, const TMatrixD& invCov, const TMatrixD& cov, int nBins);
// Scale Covariance Matrix Shape+Rate
void scaleCovShapeRate(TMatrixD& cov, const TMatrixD& covFracCollapsed, const std::vector<double>& obs, const std::vector<double>& exp, int nBins);
// Scale Covariance Matrix Shape-Only
void scaleCovShapeOnly(TMatrixD& cov, const TMatrixD& covFracCollapsed, const std::vector<double>& obs, std::vector<double>& exp, int nBins);
// Dump the Coord File:
void dumpCoordFile(int dm2_grdpts, double dm2_lowbound, double dm2_hibound, int sin22th_grdpts, double sin22th_lowbound, bool isDecay=false);
// Acquire the precomputed spectra:
void getPrecompSpec(std::vector<sbn::SBNspec>& vecSpec, int dm2_grdpts, double dm2_hibound, double dm2_lowbound, std::string xml, double scaleFactor, bool islogname = true);
// Acquire the precomputed spectra for a decay:
void getPrecompSpecDecay(std::vector<std::vector<sbn::SBNspec>>& vecSpec, int dm2_grdpts, double dm2_hibound, double dm2_lowbound, int Um4sq_grdpts, double Um4sq_hibound, double Um4sq_lowbound, std::string xml, double scaleFactor);
// Subtract Disappeared Spectra
void scaleDisappearedSpec(sbn::SBNspec& cvSpec, sbn::SBNspec massSpec, double disFactor, bool doPrint=false);
// Generate Fractional Covariance Matrices for all grid points
std::vector<std::vector<TMatrixD>> genCovMats(const sbn::SBNspec& cvSpec,
																							const std::vector<sbn::SBNspec>& a_sinsqSpec,
																							const TMatrixD& covFracSys_base,
																							const int dm2_grdpts,
																							const int sin22th_grdpts,
																							const double sin22th_lowbound,
																							const int nBins
																							);
// Create Spectra overlaid plots
void makeSpectraPlots(const sbn::SBNspec& cv, const sbn::SBNspec& app, const sbn::SBNspec& tot);
// Create Spectra Comparisons
void makeSpectraPlotsComp(const std::vector<double>& cv, const std::vector<double>& disap, int sin_base, int mi_base, double sin_val, double mi_val, int nBins=19, std::string = "exp_spectras");
void makeSpectraPlotsComp_Ratio(const std::vector<double>& cv, const std::vector<double>& disap, int sin_base, int mi_base, double sin_val, double mi_val, int nBins=19, std::string = "exp_spectras");

// Create spectra comparisons with fakedatasets and cv
void makeSpectraPlotsComp_Data(const std::vector<double>& cv, const std::vector<double>& disap, const std::vector<double>& best, int set, const std::vector<double>& unc);

// Create Exact Oscillated Spectra
sbn::SBNspec getAppSpec_exact(float mnu, bool doGen=false);



// Create nbin x nbin Empty Matrix of 0s
TMatrixD getZeroedMat(int nBins);

// }

// dM2	sin2	midval	midval
// 0.01	0.01	0.012022645	0.010964762
// 0.0144544	0.0120226	0.017378026	0.013182544
// 0.020893	0.0144544	0.02511888	0.015848929
// 0.0301995	0.017378	0.036307802	0.01905462
// 0.0436516	0.020893	0.052480742	0.022908714
// 0.0630957	0.0251189	0.075857743	0.027542299
// 0.0912011	0.0301995	0.109647965	0.0331131
// 0.131826	0.0363078	0.158489485	0.039810722
// 0.190546	0.0436516	0.229086776	0.047862997
// 0.275423	0.0524807	0.331131128	0.057543953
// 0.398107	0.0630957	0.478630016	0.069183098
// 0.57544	0.0758578	0.691831104	0.083176408
// 0.831764	0.0912011	0.999998293	0.100000091
// 1.20226	0.109648	1.44543676	0.120226691
// 1.7378	0.131826	2.089297117	0.14454401
// 2.51189	0.158489	3.019953638	0.173779875
// 3.63078	0.190546	4.365156079	0.208929681
// 5.24807	0.229087	6.309572445	0.251188831
// 7.58578	0.275423	9.12011845	0.301995188
// 10.9648	0.331131	13.18256495	0.36307791
// 15.8489	0.398107	19.05459775	0.436515697
// 22.9087	0.47863	27.54229609	0.524807438
// 33.1131	0.57544	39.81070591	0.630957392
// 47.863	0.691831	57.54398939	0.758577695
// 69.1831	0.831764	83.17637886	0.912010965
// 100	1


#endif
