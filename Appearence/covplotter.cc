#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <fstream>

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
#include "TH1F.h"
#include "TGraph.h"
// #include "TMinuitMinimizer.h"
#include "TMinuit.h"

int main(int argc, char* argv[]){
  // read in from the huge text file
  std::ifstream inFile;
	inFile.open("cov_weights.txt");
  // initialize root histogram
  TH1F* weights_h = new TH1F("weights","weights",100,0,100.000000001);
  // save weights
  float x;
  int sum =0;
  while (inFile >> x) {
    if(x>100) x=100;
    if(x>5) sum+=1;
    weights_h->Fill(x);
  }
  std::cout<<sum<<std::endl;

  TCanvas can("can", "histograms ", 1500, 1500);
  can.cd();
  weights_h->Draw("");
  can.SetLogy();
  can.SaveAs("weights.png");
  inFile.close();
	return 0;
} // end of main function
