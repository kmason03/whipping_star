#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <ctime>
#include <chrono>
// #include <gperftools/profiler.h>

#include "TFile.h"
#include "SBNllminimizer.h"
#include "prob.h"

using namespace sbn;

// define helper functions
TMatrixD GetTotalCov(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Mfracsys);
float GetLLHFromVector(const std::vector<float>& obsSpec, const SBNspec& expSpec,const TMatrixD& Msys,bool prints);
std::string ZeroPadNumber(int num, int digits=3);
std::vector<float> SetFakeData(std::vector<float> fakeData);

// define some global variables
const int nBins_e(12),nBins_mu(19);
const int nBins = nBins_e+nBins_mu;

// this function runs the FC method for a single grid point
// throws specified number of universes and returns the minimum -2LLH via a
// minimizer+ -2LLH at the point the universe was thrown from
// output is a txt file with two lines for each universe

int main(int argc, char* argv[]){

  // set the input parameter
  int specific_entry = atoi(argv[1]);

  // --------------------------------------------------
  // initiate other parameters
  // hardcoded - change to your version of the xml
  std::string sbnfithome = "/cluster/tufts/wongjiradlabnu/kmason03/taritreetests/whipping_star";
  std::string xml = sbnfithome+"/xml/TotalThreePlusOne_full.xml";
  // bins of spectra
  const int nBins_e(12),nBins_mu(19);
  const int nBins = nBins_e+nBins_mu;
  // set grid points
  const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
  const double dm2_lowbound(0.01), dm2_hibound(100);
  const double ue4_lowbound(0.01), ue4_hibound(0.5);
  const double umu4_lowbound(0.01), umu4_hibound(0.5);
  // tag for some save files
  std::string tag = "DL_full";
  // Number of universes!
  const int nFakeExp(1000);
  // --------------------------------------------------

  // initialize the minimizer
  SBNllminimizer minimizer( xml );
  std::cout<<"Initialized minimizer"<<std::endl;

  // make array to  get the parameter ids and match to the entry we want
  std::vector<std::vector<float>> params_v;
  for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
    for(int uei_base = 0; uei_base < ue4_grdpts; uei_base++){
      for(int umui_base = 0; umui_base < umu4_grdpts; umui_base++){
      	std::vector<float> params {mi_base, uei_base, umui_base};
      	params_v.push_back(params);
      }
    }
  }
  // now get the parameters we are testing
  // specific_entry=0;
  int mi_base = params_v[specific_entry][0];
  int uei_base = params_v[specific_entry][1];
  int umui_base = params_v[specific_entry][2];

  // open output files
  std::ofstream chifile;
  std::string textid = ZeroPadNumber(specific_entry, 5);
  chifile.open("chis_seek_"+textid+".txt", std::ios_base::app);


  // load in saved cvSpec - change to your location
  SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_Bkg.SBNspec.root",xml);
  // Load up the necesary bits and store them in vectors on the stack **
  cvSpec.Scale("fullosc",0.0);
  // remove mcerror, since we include it in the systematics
  cvSpec.RemoveMCError();
  std::vector<double> cv_v = cvSpec.collapsed_vector;
  // load in the full systematic matrix - created on fermilab grid
  TFile * fsys = new TFile( "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");
  TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");

  // create sbnchi - used to get cov matrix at given osc + throw universes
  SBNchi cvChi(cvSpec, *covFracSys);
  // save cv covar for plotting purposes
  TMatrixD * cvFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
  int a = cvChi.FillCollapsedFractionalMatrix(cvFracSys_collapsed);

  // there were 400 mass spectra pregenerated in the original version
  //  - this makes sure the grid points line up
  int mi_base_new = mi_base*(400/dm2_grdpts);
  //mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2
  float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
  float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
  float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
  // calculate scaling factors (sin^2(2theta) value)
  float e_app = 4*pow(ue_base,2)*pow(um_base,2);
  float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
  float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));

  // get the current test model
  std::cout << "NU Base Model: m41^2:" <<pow(mnu_base,2) <<" m41:"<<mnu_base<< " ue:" << ue_base << " um:" << um_base << std::endl;
  std::cout<<"scale factors: "<<e_app<<" "<<e_dis<<" "<<m_dis<<std::endl;
  NeutrinoModel start_model( mnu_base, uei_base, umui_base );
  start_model.Printall();
  std::cout << "mass tag: " << start_model.mass_tag << std::endl;

  minimizer._gen.regenerate_osc( start_model );
  std::cout << "Full Vector of pre-scale CV spectrum ===========" << std::endl;
  minimizer._gen.spec_central_value.PrintFullVector(true);
  std::cout << "=================================================" << std::endl;
  std::cout << "Full Vector of pre-scale dm2 spectrum ===========" << std::endl;
  minimizer._gen.spec_osc_sinsq.PrintFullVector(true);
  std::cout << "=================================================" << std::endl;

  // get the oscillated spectra as the expectation for this grid point
  SBNspec oscSpec =  minimizer.getOscSpectra( mnu_base, ue_base, um_base );
  oscSpec.RemoveMCError();

  // I'll be generating a new universe around the oscillated spectrum
  // make the chi object and fill cov matrix for this point
  SBNchi TrueChi(oscSpec, *covFracSys);
  TrueChi.ReloadCoreSpectrum(&oscSpec);
  fsys->cd();
  TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
  int b = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);

  // print statements to confirm starting spectra
  std::cout << "CV-SPEC: PrintCollapsedVector =======================" << std::endl;
  cvSpec.PrintFullVector();
  std::cout << "osc-SPEC: PrintFullVector =======================" << std::endl;
  oscSpec.PrintFullVector(true);

  // get the motor running with initial cholosky decomposition
  TrueChi.pseudo_from_collapsed = true;
  TrueChi.GeneratePseudoExperiment();

  // Across several fake experiments for this grid point:
  for(int expi = 0; expi < nFakeExp; expi++){
    std::cout<<"EXPERIMENT NUMBER: "<<expi<<std::endl;
    std::vector<float> fakeData = TrueChi.GeneratePseudoExperiment();

    // option to hardcode fakedata spectrum for testing
    if(false){
      fakeData=SetFakeData(fakeData);
    }

    // get the -2llh at the throw point (pt)
		TMatrixD cov_pt =GetTotalCov(fakeData, oscSpec, *oscFracSys_collapsed);
    float ptfit =GetLLHFromVector(fakeData, oscSpec, cov_pt, false);

    // option to run traditional grid search
    float chi_min_grid =1000000;
    if(false){
      float e_app_in = 4*pow(ue_base,2)*pow(um_base,2);
      float e_dis_in = 4*pow(ue_base,2)*(1-pow(ue_base,2));
      float m_dis_in = 4*pow(um_base,2)*(1-pow(um_base,2));
      int numgridpts =25;

      // loop over all the grid points
      for(int mi_in = 0; mi_in <numgridpts; mi_in++){
        std::cout<<mi_in<<std::endl;
        for(int uei_in = 0; uei_in < numgridpts; uei_in++){
          for(int umui_in = 0; umui_in < numgridpts; umui_in++){
            int mi_in_new = mi_in*(400/numgridpts);
            float ue_val = pow(10.,(uei_in/float(numgridpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
            float um_val = pow(10.,(umui_in/float(numgridpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
            float mnu_val = pow(10.,((mi_in+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

            // get the oscillated spectra and covariance matrix
            SBNspec inSpec =  minimizer.getOscSpectra( mnu_val, ue_val, um_val );
            // inSpec.PrintCollapsedVector();
            inSpec.RemoveMCError();
            SBNchi innerChi(inSpec, *covFracSys);
            TMatrixD * inFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
            int c = innerChi.FillCollapsedFractionalMatrix(inFracSys_collapsed);
            TMatrixD cov_grid = GetTotalCov(fakeData, inSpec, *inFracSys_collapsed);
            float chi_test = GetLLHFromVector(fakeData, inSpec, cov_grid, false);
            delete inFracSys_collapsed;

            if (chi_test<chi_min_grid){
              chi_min_grid = chi_test;
            }
          }
        }
      }
    }

    // now run the minimizer
    double bestfit = minimizer.doFit( fakeData, mnu_base, ue_base, um_base )[0];
    // double bestfit = 0;


    // trying a variety of starts
    double min_minimizer =100000000;
    std::vector<float> beststart;
    if(true){
      std::cout<<"starting loop"<<std::endl;
      double temp;
      // std::vector<float> m_v={0.01,2,8};
      // std::vector<float> u_v={0.01,0.08,0.4};
      for(int x=0;x<5;x++){
            int startval = 5*x;
            int numgridpts=25;
            float ue_val = pow(10.,(startval/float(numgridpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
            float um_val = pow(10.,(startval/float(numgridpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
            float mnu_val = pow(10.,((startval*(500/float(numgridpts))+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

            temp = minimizer.doFit( fakeData, mnu_val,ue_val,um_val)[0];
            if (temp < min_minimizer){
              min_minimizer=temp;
              // beststart={m_v[x],u_v[y],u_v[z]};
          //   }
          // }
        }
      } //end of loop over start points
    }

    chifile<<ptfit<<std::endl;
    // chifile<<bestfit<<std::endl;
    chifile<<min_minimizer<<std::endl;
    // chifile<<chi_min_grid<<std::endl;
  } //end of fake experiment loop
  return 0;
} // end of main function

TMatrixD GetTotalCov(const std::vector<float>& obsSpec, const SBNspec& expSpec, const TMatrixD& Mfracsys){
	// function to take the fractional Msys and return totl Msys+Mstat
	// inputs:
	// obsSpec: "data" spectra
	// predSpec: "MC" spectra
	// Mfracsys: fractional (flux+xsec+detvar) covariance matrix
	TMatrixD fullcov(nBins,nBins);
	for(int i = 0; i < nBins; i++){
		for(int j = 0; j < nBins; j++){
			// first set to zero
			fullcov[i][j] = 0.0;
			// scale to the prediction
			fullcov[i][j] = (Mfracsys)[i][j]*expSpec.collapsed_vector[i]*expSpec.collapsed_vector[j];
			// add in stat errors start with CNP for "data" errors;
			if(i==j){
				if (expSpec.collapsed_vector[i] >0){
					// 1/observed 2/expectation
					fullcov[i][j] += 3.0 / (1.0/obsSpec[i] + 2.0/expSpec.collapsed_vector[i]);
				}
				else {
					fullcov[i][j] += expSpec.collapsed_vector[i]/2.0;
				}
			}
		}//end of first bin loop
	}//end of second bin loop
	return fullcov;
}//end of GetTotalCov


float GetLLHFromVector(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Msys, bool prints){
	// // function to calculate a chi2 (shape + rate)
	// // inputs:
	// // obsSpec: "data" vector
	// // predSpec: "MC" spectra
	// // Mfracsys: total (flux+xsec+detvar) covariance matrix
  float chisqTest = 0;
  // inv cov
  TMatrixD invcov = Msys;
  invcov.Invert();

  // add the chi2-like part
  chisqTest = 0;
  for(int i = 0; i < nBins; i++){
    for(int j = 0; j < nBins; j++){
    if(i==j && false) std::cout<<i<<" "
              <<obsSpec[i]<<" "
              <<expSpec.collapsed_vector[i]<<" "
              <<Msys[i][j]<<" "
              <<((obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]))
              <<std::endl;
    chisqTest += (obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]);
        }
      }
  // now need ln(det(2Pi*M))
  chisqTest += log(Msys.Determinant());

  return chisqTest;
}//end of GetLLHFromSpectra

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

std::vector<float> SetFakeData(std::vector<float> fakeData){

  // these are the data values!
  fakeData[0] = 4;
  fakeData[1] = 1;
  fakeData[2] = 1;
  fakeData[3] = 2;
  fakeData[4] = 5;
  fakeData[5] = 3;
  fakeData[6] = 8;
  fakeData[7] = 0;
  fakeData[8] = 1;
  fakeData[9] = 0;
  fakeData[10] = 4;
  fakeData[11] = 2;
  fakeData[12] = 26;
  fakeData[13] = 192;
  fakeData[14] = 276;
  fakeData[15] = 401;
  fakeData[16] = 389;
  fakeData[17] = 463;
  fakeData[18] = 439;
  fakeData[19] = 482;
  fakeData[20] = 395;
  fakeData[21] = 353;
  fakeData[22] = 303;
  fakeData[23] = 233;
  fakeData[24] = 240;
  fakeData[25] = 177;
  fakeData[26] = 118;
  fakeData[27] = 109;
  fakeData[28] = 109;
  fakeData[29] = 85;
  fakeData[30] = 58;
  return fakeData;
}
