#include <array>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <boost/algorithm/string.hpp>

int main() {

  // Initialize names of input files
  // ... for fake data event lists
  const char* sel_run1_eventlist_file = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_data_run1_eventlist.txt";
  const char* sel_run2_eventlist_file = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_data_run2_eventlist.txt";
  const char* sel_run3_eventlist_file = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_data_run3_eventlist.txt";
  // ... for output
  const char* output_file = "/uboone/data/users/yatesla/othersys_mcc9/input_to_sbnfit/input_to_sbnfit_data_1mu1p_Jun28.root";
  // ... for flagging whether there is a vtxid between the RSE and reco energy
  const bool has_vtxid = true;

  // Initalize the number of reco variables
  const int N_var = 75;
  
  //** Read in lists of selected events and their reconstructed energies **//
  
  // Initialize maps to hold selected events
  std::map<std::array<int,3>,std::array<double,N_var>> sel_run1_event_map;
  std::map<std::array<int,3>,std::array<double,N_var>> sel_run2_event_map;
  std::map<std::array<int,3>,std::array<double,N_var>> sel_run3_event_map;
  
  // Initialize variables to hold intermediate values
  std::string line;
  std::vector<std::string> split_line;
  std::array<int,3> rse;
  std::array<double,N_var> reco_var;

  // Read in the Run 1 event list
  std::ifstream sel_run1_eventlist(sel_run1_eventlist_file);
  if ( sel_run1_eventlist.is_open() ) {
    while ( getline(sel_run1_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      //std::cout << rse[0] << ", " << rse[1] << ", " << rse[2] << std::endl;
      for ( int i=0; i<N_var; i++ ) {
	//std:: cout << i << ": ";
	if ( has_vtxid ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
	else reco_var[i] = std::stof(split_line[i+3]);
	//std::cout << reco_var[i] << std::endl;
      }
      sel_run1_event_map.emplace(rse,reco_var);
    }
    sel_run1_eventlist.close();
  }
  std::cout << "Read in " << sel_run1_event_map.size() << " selected events from the Run 1 sample" << std::endl;

  // Read in the Run 2 event list
  std::ifstream sel_run2_eventlist(sel_run2_eventlist_file);
  if ( sel_run2_eventlist.is_open() ) {
    while ( getline(sel_run2_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<N_var; i++ ) {
	if ( has_vtxid ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
	else reco_var[i] = std::stof(split_line[i+3]);
      }
      sel_run2_event_map.emplace(rse,reco_var);
    }
    sel_run2_eventlist.close();
  }
  std::cout << "Read in " << sel_run2_event_map.size() << " selected events from the Run 2 sample" << std::endl;

  // Read in the Run 3 event list
  std::ifstream sel_run3_eventlist(sel_run3_eventlist_file);
  if ( sel_run3_eventlist.is_open() ) {
    while ( getline(sel_run3_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<N_var; i++ ) {
	if ( has_vtxid ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
	else reco_var[i] = std::stof(split_line[i+3]);
      }
      sel_run3_event_map.emplace(rse,reco_var);
    }
    sel_run3_eventlist.close();
  }
  std::cout << "Read in " << sel_run3_event_map.size() << " selected events from the Run 3 sample" << std::endl;
  
  
  // * Initialize trees in the the output file * //
  
  // Open the output file
  TFile* out = new TFile(output_file, "RECREATE");
  if ( !out->IsOpen() ) std::cout << "Failed to open file " << output_file << std::endl;
  
  // Initialize variables for the desired branches
  int run, subrun, event;
  //double nu_energy_reco;  // defined as the first entry in reco_var, initialized above
  std::map<std::string, std::vector<double>> sys_weights;
  
  // Initialize output trees
  // ... for Run 1
  TTree* sel_run1_output_tree = new TTree("sel_run1_tree", "events selected from Run 1 sample");
  sel_run1_output_tree->Branch("run",&run);
  sel_run1_output_tree->Branch("subrun",&subrun);
  sel_run1_output_tree->Branch("event",&event);
  sel_run1_output_tree->Branch("weights",&sys_weights);
  sel_run1_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_run1_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_run1_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_run1_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_run1_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_run1_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_run1_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_run1_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_run1_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_run1_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_run1_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_run1_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_run1_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_run1_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_run1_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_run1_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_run1_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_run1_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_run1_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_run1_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_run1_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_run1_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_run1_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_run1_output_tree->Branch("x_reco",&reco_var[23]);
  sel_run1_output_tree->Branch("y_reco",&reco_var[24]);
  sel_run1_output_tree->Branch("z_reco",&reco_var[25]);
  sel_run1_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_run1_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_run1_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_run1_output_tree->Branch("phiT_reco",&reco_var[29]);
  sel_run1_output_tree->Branch("Q2_reco",&reco_var[30]);
  sel_run1_output_tree->Branch("lepton_length_reco",&reco_var[31]);
  sel_run1_output_tree->Branch("proton_length_reco",&reco_var[32]);
  sel_run1_output_tree->Branch("proton_cos_theta_reco",&reco_var[33]);
  sel_run1_output_tree->Branch("lepton_cos_theta_reco",&reco_var[34]);
  // skip reco_var[35-37], which are dummy values for true variables
  sel_run1_output_tree->Branch("dllee_pass_simple_cuts_reco",&reco_var[38]);
  sel_run1_output_tree->Branch("dllee_failed_boost_reco",&reco_var[39]);
  sel_run1_output_tree->Branch("tot_PE_reco",&reco_var[40]);
  sel_run1_output_tree->Branch("porch_tot_PE_reco",&reco_var[41]);
  // skip reco_var[42-46], which are also dummy values
  sel_run1_output_tree->Branch("dllee_nBDTs",&reco_var[47]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[48]);
  sel_run1_output_tree->Branch("dllee_bdt_score00",&reco_var[49]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[50]);
  sel_run1_output_tree->Branch("dllee_bdt_score01",&reco_var[51]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[52]);
  sel_run1_output_tree->Branch("dllee_bdt_score02",&reco_var[53]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[54]);
  sel_run1_output_tree->Branch("dllee_bdt_score03",&reco_var[55]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[56]);
  sel_run1_output_tree->Branch("dllee_bdt_score04",&reco_var[57]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[58]);
  sel_run1_output_tree->Branch("dllee_bdt_score05",&reco_var[59]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[60]);
  sel_run1_output_tree->Branch("dllee_bdt_score06",&reco_var[61]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[62]);
  sel_run1_output_tree->Branch("dllee_bdt_score07",&reco_var[63]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[64]);
  sel_run1_output_tree->Branch("dllee_bdt_score08",&reco_var[65]);
  sel_run1_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[66]);
  sel_run1_output_tree->Branch("dllee_bdt_score09",&reco_var[67]);
  sel_run1_output_tree->Branch("dllee_bdt_score_avg",&reco_var[68]);
  sel_run1_output_tree->Branch("dllee_bdt_score_median",&reco_var[69]);
  sel_run1_output_tree->Branch("dllee_bdt_score_max",&reco_var[70]);
  // skip reco_var[71] and reco_var[72], which are not clear to me
  sel_run1_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[73]);
  sel_run1_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[74]);
  // ... for Run 2
  TTree* sel_run2_output_tree = new TTree("sel_run2_tree", "events selected from Run 2 sample");
  sel_run2_output_tree->Branch("run",&run);
  sel_run2_output_tree->Branch("subrun",&subrun);
  sel_run2_output_tree->Branch("event",&event);
  sel_run2_output_tree->Branch("weights",&sys_weights);
  sel_run2_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_run2_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_run2_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_run2_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_run2_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_run2_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_run2_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_run2_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_run2_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_run2_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_run2_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_run2_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_run2_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_run2_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_run2_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_run2_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_run2_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_run2_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_run2_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_run2_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_run2_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_run2_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_run2_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_run2_output_tree->Branch("x_reco",&reco_var[23]);
  sel_run2_output_tree->Branch("y_reco",&reco_var[24]);
  sel_run2_output_tree->Branch("z_reco",&reco_var[25]);
  sel_run2_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_run2_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_run2_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_run2_output_tree->Branch("phiT_reco",&reco_var[29]);
  sel_run2_output_tree->Branch("Q2_reco",&reco_var[30]);
  sel_run2_output_tree->Branch("lepton_length_reco",&reco_var[31]);
  sel_run2_output_tree->Branch("proton_length_reco",&reco_var[32]);
  sel_run2_output_tree->Branch("proton_cos_theta_reco",&reco_var[33]);
  sel_run2_output_tree->Branch("lepton_cos_theta_reco",&reco_var[34]);
  sel_run2_output_tree->Branch("dllee_pass_simple_cuts_reco",&reco_var[38]);
  sel_run2_output_tree->Branch("dllee_failed_boost_reco",&reco_var[39]);
  sel_run2_output_tree->Branch("tot_PE_reco",&reco_var[40]);
  sel_run2_output_tree->Branch("porch_tot_PE_reco",&reco_var[41]);
  sel_run2_output_tree->Branch("dllee_nBDTs",&reco_var[47]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[48]);
  sel_run2_output_tree->Branch("dllee_bdt_score00",&reco_var[49]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[50]);
  sel_run2_output_tree->Branch("dllee_bdt_score01",&reco_var[51]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[52]);
  sel_run2_output_tree->Branch("dllee_bdt_score02",&reco_var[53]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[54]);
  sel_run2_output_tree->Branch("dllee_bdt_score03",&reco_var[55]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[56]);
  sel_run2_output_tree->Branch("dllee_bdt_score04",&reco_var[57]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[58]);
  sel_run2_output_tree->Branch("dllee_bdt_score05",&reco_var[59]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[60]);
  sel_run2_output_tree->Branch("dllee_bdt_score06",&reco_var[61]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[62]);
  sel_run2_output_tree->Branch("dllee_bdt_score07",&reco_var[63]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[64]);
  sel_run2_output_tree->Branch("dllee_bdt_score08",&reco_var[65]);
  sel_run2_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[66]);
  sel_run2_output_tree->Branch("dllee_bdt_score09",&reco_var[67]);
  sel_run2_output_tree->Branch("dllee_bdt_score_avg",&reco_var[68]);
  sel_run2_output_tree->Branch("dllee_bdt_score_median",&reco_var[69]);
  sel_run2_output_tree->Branch("dllee_bdt_score_max",&reco_var[70]);
  sel_run2_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[73]);
  sel_run2_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[74]);
  // ... for Run 3
  TTree* sel_run3_output_tree = new TTree("sel_run3_tree", "events selected from Run 3 sample");
  sel_run3_output_tree->Branch("run",&run);
  sel_run3_output_tree->Branch("subrun",&subrun);
  sel_run3_output_tree->Branch("event",&event);
  sel_run3_output_tree->Branch("weights",&sys_weights);
  sel_run3_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_run3_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_run3_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_run3_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_run3_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_run3_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_run3_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_run3_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_run3_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_run3_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_run3_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_run3_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_run3_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_run3_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_run3_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_run3_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_run3_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_run3_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_run3_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_run3_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_run3_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_run3_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_run3_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_run3_output_tree->Branch("x_reco",&reco_var[23]);
  sel_run3_output_tree->Branch("y_reco",&reco_var[24]);
  sel_run3_output_tree->Branch("z_reco",&reco_var[25]);
  sel_run3_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_run3_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_run3_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_run3_output_tree->Branch("phiT_reco",&reco_var[29]);
  sel_run3_output_tree->Branch("Q2_reco",&reco_var[30]);
  sel_run3_output_tree->Branch("lepton_length_reco",&reco_var[31]);
  sel_run3_output_tree->Branch("proton_length_reco",&reco_var[32]);
  sel_run3_output_tree->Branch("proton_cos_theta_reco",&reco_var[33]);
  sel_run3_output_tree->Branch("lepton_cos_theta_reco",&reco_var[34]);
  sel_run3_output_tree->Branch("dllee_pass_simple_cuts_reco",&reco_var[38]);
  sel_run3_output_tree->Branch("dllee_failed_boost_reco",&reco_var[39]);
  sel_run3_output_tree->Branch("tot_PE_reco",&reco_var[40]);
  sel_run3_output_tree->Branch("porch_tot_PE_reco",&reco_var[41]);
  sel_run3_output_tree->Branch("dllee_nBDTs",&reco_var[47]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[48]);
  sel_run3_output_tree->Branch("dllee_bdt_score00",&reco_var[49]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[50]);
  sel_run3_output_tree->Branch("dllee_bdt_score01",&reco_var[51]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[52]);
  sel_run3_output_tree->Branch("dllee_bdt_score02",&reco_var[53]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[54]);
  sel_run3_output_tree->Branch("dllee_bdt_score03",&reco_var[55]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[56]);
  sel_run3_output_tree->Branch("dllee_bdt_score04",&reco_var[57]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[58]);
  sel_run3_output_tree->Branch("dllee_bdt_score05",&reco_var[59]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[60]);
  sel_run3_output_tree->Branch("dllee_bdt_score06",&reco_var[61]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[62]);
  sel_run3_output_tree->Branch("dllee_bdt_score07",&reco_var[63]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[64]);
  sel_run3_output_tree->Branch("dllee_bdt_score08",&reco_var[65]);
  sel_run3_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[66]);
  sel_run3_output_tree->Branch("dllee_bdt_score09",&reco_var[67]);
  sel_run3_output_tree->Branch("dllee_bdt_score_avg",&reco_var[68]);
  sel_run3_output_tree->Branch("dllee_bdt_score_median",&reco_var[69]);
  sel_run3_output_tree->Branch("dllee_bdt_score_max",&reco_var[70]);
  sel_run3_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[73]);
  sel_run3_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[74]);
  
  std::cout << "Initialized trees in output files" << std::endl;
  
  
  // * Fill trees in the output file * //

  // Fill Run 1 tree
  for ( auto event_it = sel_run1_event_map.begin(); event_it != sel_run1_event_map.end(); event_it++ ) {
    // Get the variables to be put into the tree                                                                                                                                 
    run = event_it->first[0];
    subrun = event_it->first[1];
    event = event_it->first[2];
    reco_var = event_it->second;
    // Fill the tree                                                                                                                                  
    sel_run1_output_tree->Fill();
  }

  // Fill Run 2 tree
  for ( auto event_it = sel_run2_event_map.begin(); event_it != sel_run2_event_map.end(); event_it++ ) {
    // Get the variables to be put into the tree                                                                                                                                 
    run = event_it->first[0];
    subrun = event_it->first[1];
    event = event_it->first[2];
    reco_var = event_it->second;
    // Fill the tree                                                                                                                                  
    sel_run2_output_tree->Fill();
  }
  
  // Fill Run 3 tree
  for ( auto event_it = sel_run3_event_map.begin(); event_it != sel_run3_event_map.end(); event_it++ ) {
    // Get the variables to be put into the tree                                                                                                                                 
    run = event_it->first[0];
    subrun = event_it->first[1];
    event = event_it->first[2];
    reco_var = event_it->second;
    // Fill the tree                                                                                                                                                             
    sel_run3_output_tree->Fill();
  }

  std::cout << "Filled output trees" << std::endl;


  // * Conclusion * //
  
  // Count events in all output trees and make sure we did this right
  int n_sel_run1_in  = sel_run1_event_map.size();
  int n_sel_run1_out = sel_run1_output_tree->GetEntries();
  int n_sel_run2_in  = sel_run2_event_map.size();
  int n_sel_run2_out = sel_run2_output_tree->GetEntries();
  int n_sel_run3_in  = sel_run3_event_map.size();
  int n_sel_run3_out = sel_run3_output_tree->GetEntries();
  
  // Print information
  std::cout << "Selected event count information" << std::endl;
  std::cout << "  Run 1 events: " << n_sel_run1_in << " events in, " << n_sel_run1_out << " events out" << std::endl;
  std::cout << "  Run 2 events: " << n_sel_run2_in << " events in, " << n_sel_run2_out << " events out" << std::endl;
  std::cout << "  Run 3 events: " << n_sel_run3_in << " events in, " << n_sel_run3_out << " events out" << std::endl;
  if ( n_sel_run1_in==n_sel_run1_out and  n_sel_run2_in==n_sel_run2_out and n_sel_run3_in==n_sel_run3_out ) {
    std::cout << "These values are consistent!" << std::endl;
  }
  else {
    std::cout << "Warning: These values are NOT consistent!" << std::endl;
  }
  
  // Write and close the output file
  out->Write();
  out->Close();
  
  // El Fin
  std::cout << "Wrote output file " << output_file << std::endl;
  std::cout << "Done!" << std::endl;
  return 0;
  
}
