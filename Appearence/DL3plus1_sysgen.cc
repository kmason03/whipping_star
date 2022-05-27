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
void generate_sys(int massid);

// define some global variables
std::string xml = "/uboone/app/users/kmason/whipping_star/xml/TotalThreePlusOne_full.xml";
std::string tag = "DL_full";
// set these parameters at the very start
const double dm2_lowbound(0.01), dm2_hibound(100);
const int dm2_grdpts(25);
double  mnu;


int main(int argc, char* argv[]){
  generate_sys();
  return 0;

} // end of main function


//-----------------------------------------------------------------------------
//----------------------HELPER FUNCTIONS---------------------------------------
//----------------------------------------------------------------------------

 void generate_sys(){
	 //function to generate the flux+xsec+reweight systematic covariance matrix
	 // prerunning this speeds up the initial grid search part and helps with plotting
	 // note: only works on fermilab
	 // writes covar files to current directory

	// start with a null model if mass id ==0
     SBNcovariance _covar(xml);
    _covar.FormCovarianceMatrix(tag);
    //_covar.PrintVariations(tag);
    _covar.PrintMatricies(tag);
    _covar.frac_covariance.Print();
    NeutrinoModel nullModel(0, 0, 0);
    SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
    SBNspec bkg = bkgo->spec_central_value;
    bkg.WriteOut(tag+"_Bkg");

    return;
}//end of generate spectra function
