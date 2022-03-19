#include <iostream>
#include <string>
#include <chrono>
#include "SBNgenerate.h"
#include "prob.h"


using namespace sbn;

// define some global variables

int main(int argc, char* argv[]){

  // NeutrinoModel nullModel(0, 0, 0);
  // SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
  // SBNspec bkg = bkgo->spec_central_value;
  
  std::string xml = "../xml/TotalThreePlusOne_full.xml";
  std::string tag = "DL_full";
  // set these parameters at the very start
  const double dm2_lowbound(0.01), dm2_hibound(100);
  const int dm2_grdpts(400);
  double  mnu;
  
  //we're precomputing, so we don't really care about the u's set to sin2=1
  float dm2_eV2 = 2.0; // eV2
  float dm_eV = sqrt( dm2_eV2 );
  
  //Model: mnu, ue4, um4
  //mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
  //see davio's thesis eqs 4.14-16
  
  NeutrinoModel testModel(dm_eV, sqrt(.5), sqrt(.5));

  // on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
  SBNgenerate * gen = new SBNgenerate(xml,testModel, true);

  // save the initial spectrum
  auto const& branch_variable = gen->branch_variables[0][0];
  int ih = gen->spec_central_value.map_hist.at(branch_variable->associated_hist);
  std::cout << "hist ih=" << ih << std::endl;
  std::map< std::string, TH1D* > hists;
  for ( auto it=gen->spec_central_value.map_hist.begin(); it!=gen->spec_central_value.map_hist.end(); it++ ) {
    std::string name = it->first;
    std::string name_copy = name+"_clone";
    TH1D* hist_copy = (TH1D*) gen->spec_central_value.hist.at(it->second).Clone( name_copy.c_str() );
    hists[ name ] = hist_copy;
  }

  // for(int j=0;j<gen->num_files;j++){
  //   std::cout << "FILE[" << j << "] ----------------------------" << std::endl;
  //   std::cout << "  num branch vars: " << gen->branch_variables[j].size() << std::endl;
  //   for(int t=0; t<gen->branch_variables[j].size();t++){
  //     auto const& b_var = gen->branch_variables[j][t];      
  //     ih = gen->spec_central_value.map_hist.at(b_var->associated_hist);
  //     std::cout << "[" << j << "][" << t << "] hist=" << gen->spec_central_value.hist[ ih ].GetName()
  // 		<< " Nbins=" << gen->spec_central_value.hist[ ih ].GetXaxis()->GetNbins()
  // 		<< std::endl;
  //   }
  // }

  for ( auto it=hists.begin(); it!=hists.end(); it++ ) {
    std::cout << it->first << "-----------------" << std::endl;
    for (int ib=1; ib<=it->second->GetXaxis()->GetNbins(); ib++) {
      std::cout << "bin[" << ib << "] "
		<< " cv=" << it->second->GetBinContent(ib)
	//<< " osc=" << gen->spec_osc_sinsq.map_hist[it->first].GetBinContent(ib)
		<< std::endl;
    }
  }
  
  for (auto const& cache : gen->event_cache_v ) {
    std::cout << "in cache: " << cache.num_events << std::endl;
  }
  
  std::clock_t start = std::clock();

  int num_refills = 100;

  for (int i=0; i<num_refills; i++) {
    gen->regenerate_osc( testModel );
    
    if (i==0) {
      // confirm the regen reproduces previous histograms
      float totdiff = 0.;
      for ( auto it=gen->spec_central_value.map_hist.begin(); it!=gen->spec_central_value.map_hist.end(); it++ ) {
	std::string name = it->first;
	TH1D* orig = hists[name];
	auto const& h = gen->spec_central_value.hist[ it->second ];
	for (int ib=1; ib<=h.GetXaxis()->GetNbins(); ib++) {
	  totdiff += fabs( h.GetBinContent(ib) - orig->GetBinContent(ib) );
	}
      }
      std::cout << "Total Diff: " << totdiff << std::endl;
    }

  }//end of refill loop
  
  std::clock_t end   = std::clock();

  float dt = ((float)end - start)/CLOCKS_PER_SEC/(float)num_refills;
  std::cout << "time per regen: " << dt << std::endl;  
  
  return 0;

} // end of main function


