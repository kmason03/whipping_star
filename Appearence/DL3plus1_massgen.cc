#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <cstdlib>

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
// #include "TMinuitMinimizer.h"
#include "TMinuit.h"

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

// define helper functions
void generate_spectra(int massid);

// define some global variables
std::string xml = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/xml/TotalThreePlusOne_full.xml";
std::string tag = "DL_full";
// set these parameters at the very start
const double dm2_lowbound(0.01), dm2_hibound(100);
const int dm2_grdpts(400);
double  mnu;


int main(int argc, char* argv[]){
  // get the input integer: require 1
  // this script pregenerates spectra for the given delta m^2 - not up to date with code speed up,
  // used to make thesis plots
  
  int specific_entry = atoi(argv[1]);
  std::cout<<specific_entry<<std::endl;
	generate_spectra(specific_entry);
	return 0;
  // * note you only need tag_SINSQ.root files. Can delete the others when done
} // end of main function


//-----------------------------------------------------------------------------
//----------------------HELPER FUNCTIONS---------------------------------------
//----------------------------------------------------------------------------

 void generate_spectra(int massid){
	 //function to generate the different mass spectra we need
	 // prerunning this speeds up the initial grid search part and helps with plotting
	 // note: currently run twice to prevent the job being killed
	 // inputs:
	 // int massid: mass index to run
	 // outputs:
	 // none, but writes spectra in root files to current directory

	// start with a null model if mass id ==0
  if(massid==0){
    NeutrinoModel nullModel(0, 0, 0);
    SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
    SBNspec bkg = bkgo->spec_central_value;
    bkg.WriteOut(tag+"_Bkg");
  }

	mnu = pow(10.,((massid+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
	//Model: mnu, ue4, um4
	//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
	//see davio's thesis eqs 4.14-16
	//we're precomputing, so we don't really care about the u's set to sin2=1
	NeutrinoModel testModel(mnu, sqrt(.5), sqrt(.5));

	// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
	SBNgenerate * gen = new SBNgenerate(xml,testModel);
	std::cout<<massid<<std::endl;
	// Write them to files
	gen->WritePrecomputedOscSpecs(tag);

	return;
}//end of generate spectra function
