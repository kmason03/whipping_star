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
  // ... for event lists
  const char* sel1e1p_nue_eventlist_file     = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1e1p_nue_eventlist.txt";
  const char* sel1e1p_bnb_eventlist_file     = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1e1p_bnb_eventlist.txt";
  const char* sel1e1p_dirt_eventlist_file    = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1e1p_dirt_eventlist.txt";
  const char* sel1e1p_extbnb_eventlist_file  = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1e1p_extbnb_eventlist.txt";
  const char* sel1mu1p_bnb_eventlist_file    = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_bnb_eventlist.txt";
  const char* sel1mu1p_nue_eventlist_file    = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_nue_eventlist.txt";
  const char* sel1mu1p_dirt_eventlist_file   = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_dirt_eventlist.txt";
  const char* sel1mu1p_extbnb_eventlist_file = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_extbnb_eventlist.txt";
  // ... for 5e19 event list
  const char* sel1mu1p_5e19_eventlist_file = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_5e19_eventlist.txt";
  // ... for Tune 1 weights
  const char* v304_xsec_spline_file  = "/uboone/data/users/yatesla/othersys_mcc9/for_Tune1_weights/xsec_graphs_mcc9_v304.root";
  const char* tune1_xsec_spline_file = "/uboone/data/users/yatesla/othersys_mcc9/for_Tune1_weights/xsec_graphs_tune1.root";
  // ... for arborist files
  const char* nue_arborist_file  = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v40_intrinsic_nue_run1.root";
  const char* bnb_arborist_file  = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v40_bnb_nu_run1.root";
  const char* dirt_arborist_file = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v40_dirt_nu_run1.root";
  // ... for output
  const char* output_file = "/uboone/data/users/yatesla/othersys_mcc9/input_to_sbnfit/input_to_sbnfit_v40_Apr27.root";
  
  
  //** Read in lists of selected events and their reconstructed energies **//
  
  // Initialize maps to hold selected events
  std::map<std::array<int,3>,double> sel1e1p_nue_event_map;
  std::map<std::array<int,3>,double> sel1e1p_bnb_event_map;
  std::map<std::array<int,3>,double> sel1e1p_dirt_event_map;
  std::map<std::array<int,3>,double> sel1e1p_extbnb_event_map;
  std::map<std::array<int,3>,std::array<double,29>> sel1mu1p_bnb_event_map;
  std::map<std::array<int,3>,std::array<double,29>> sel1mu1p_nue_event_map;
  std::map<std::array<int,3>,std::array<double,29>> sel1mu1p_dirt_event_map;
  std::map<std::array<int,3>,std::array<double,29>> sel1mu1p_extbnb_event_map;
  std::map<std::array<int,3>,double> sel1mu1p_5e19_event_map;
  
  // Initialize variables to hold intermediate values
  std::string line;
  std::vector<std::string> split_line;
  std::array<int,3> rse;
  double reco_energy;
  std::array<double,29> reco_var;
  
  // Read in the nue intrinsic event list
  std::ifstream sel1e1p_nue_eventlist(sel1e1p_nue_eventlist_file);
  if ( sel1e1p_nue_eventlist.is_open() ) {
    while ( getline(sel1e1p_nue_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      reco_energy = std::stof(split_line[3]);
      sel1e1p_nue_event_map.emplace(rse,reco_energy);
    }
    sel1e1p_nue_eventlist.close();
  }
  std::cout << "Read in " << sel1e1p_nue_event_map.size() << " selected 1e1p events from the intinsic nue overlay sample" << std::endl;
  
  // Read in the nue bnb event list
  std::ifstream sel1e1p_bnb_eventlist(sel1e1p_bnb_eventlist_file);
  if ( sel1e1p_bnb_eventlist.is_open() ) {
    while ( getline(sel1e1p_bnb_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      reco_energy = std::stof(split_line[3]);
      sel1e1p_bnb_event_map.emplace(rse,reco_energy);
    }
    sel1e1p_bnb_eventlist.close();
  }
  std::cout << "Read in " << sel1e1p_bnb_event_map.size() << " selected 1e1p events from the BNB overlay sample" << std::endl;
  
  // Read in the nue dirt event list
  std::ifstream sel1e1p_dirt_eventlist(sel1e1p_dirt_eventlist_file);
  if ( sel1e1p_dirt_eventlist.is_open() ) {
    while ( getline(sel1e1p_dirt_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      reco_energy = std::stof(split_line[3]);
      sel1e1p_dirt_event_map.emplace(rse,reco_energy);
    }
    sel1e1p_dirt_eventlist.close();
  }
  std::cout << "Read in " << sel1e1p_dirt_event_map.size() << " selected 1e1p events from the dirt overlay sample" << std::endl;
  
  // Read in the nue extbnb event list
  std::ifstream sel1e1p_extbnb_eventlist(sel1e1p_extbnb_eventlist_file);
  if ( sel1e1p_extbnb_eventlist.is_open() ) {
    while ( getline(sel1e1p_extbnb_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      reco_energy = std::stof(split_line[3]);
      sel1e1p_extbnb_event_map.emplace(rse,reco_energy);
    }
    sel1e1p_extbnb_eventlist.close();
  }
  std::cout << "Read in " << sel1e1p_extbnb_event_map.size() << " selected 1e1p events from the EXTBNB data sample" << std::endl;
  
  // Read in numu bnb event list
  std::ifstream sel1mu1p_bnb_eventlist(sel1mu1p_bnb_eventlist_file);
  if ( sel1mu1p_bnb_eventlist.is_open() ) {
    while ( getline(sel1mu1p_bnb_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<29; i++ ) reco_var[i] = std::stof(split_line[i+3]);
      sel1mu1p_bnb_event_map.emplace(rse,reco_var);
      //std::cout << rse[0] << ", " << rse[1] << ", " << rse[2] << std::endl;
    }
    sel1mu1p_bnb_eventlist.close();
  }
  std::cout << "Read in " << sel1mu1p_bnb_event_map.size() << " selected 1mu1p events from the BNB overlay sample" << std::endl;

  // Read in numu nue event list
  std::ifstream sel1mu1p_nue_eventlist(sel1mu1p_nue_eventlist_file);
  if ( sel1mu1p_nue_eventlist.is_open() ) {
    while ( getline(sel1mu1p_nue_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<29; i++ ) reco_var[i] = std::stof(split_line[i+3]);
      sel1mu1p_nue_event_map.emplace(rse,reco_var);
    }
    sel1mu1p_nue_eventlist.close();
  }
  std::cout << "Read in " << sel1mu1p_nue_event_map.size() << " selected 1mu1p events from the intrinsic nue overlay sample" << std::endl;
  
  // Read in the numu dirt event list
  std::ifstream sel1mu1p_dirt_eventlist(sel1mu1p_dirt_eventlist_file);
  if ( sel1mu1p_dirt_eventlist.is_open() ) {
    while ( getline(sel1mu1p_dirt_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<29; i++ ) reco_var[i] = std::stof(split_line[i+3]);
      sel1mu1p_dirt_event_map.emplace(rse,reco_var);
    }
    sel1mu1p_dirt_eventlist.close();
  }
  std::cout << "Read in " << sel1mu1p_dirt_event_map.size() << " selected 1mu1p events from the dirt overlay sample" << std::endl;
  
  // Read in the numu extbnb event list
  std::ifstream sel1mu1p_extbnb_eventlist(sel1mu1p_extbnb_eventlist_file);
  if ( sel1mu1p_extbnb_eventlist.is_open() ) {
    while ( getline(sel1mu1p_extbnb_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<29; i++ ) reco_var[i] = std::stof(split_line[i+3]);
      sel1mu1p_extbnb_event_map.emplace(rse,reco_var);
    }
    sel1mu1p_extbnb_eventlist.close();
  }
  std::cout << "Read in " << sel1mu1p_extbnb_event_map.size() << " selected 1mu1p events from the EXTBNB data sample" << std::endl;
  
  // Read in the numu bnb5e19 event list
  std::ifstream sel1mu1p_5e19_eventlist(sel1mu1p_5e19_eventlist_file);
  if ( sel1mu1p_5e19_eventlist.is_open() ) {
    while ( getline(sel1mu1p_5e19_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      reco_energy = std::stof(split_line[3]);
      sel1mu1p_5e19_event_map.emplace(rse,reco_energy);
    }
    sel1mu1p_5e19_eventlist.close();
  }
  std::cout << "Read in " << sel1mu1p_5e19_event_map.size() << " selected 1mu1p events from the 5e19 data sample" << std::endl;
  
  
  // * Initialize trees in the the output file * //
  
  // Open the output file
  TFile* out = new TFile(output_file, "RECREATE");
  if ( !out->IsOpen() ) std::cout << "Failed to open file " << output_file << std::endl;
  
  // Initialize variables for the desired branches
  int run, subrun, event;
  double nu_energy_reco;
  int nu_pdg;
  double nu_energy_true, nu_L_true;
  int nu_interaction_ccnc, nu_interaction_mode, nu_interaction_type, nu_target_pdg;
  double tune1_weight, spline_weight, rootino_weight, ub_tune_weight, xsec_corr_weight, lee_weight;
  std::map<std::string, std::vector<double>> sys_weights;
  auto* sys_weights_ptr = &sys_weights;
  // reco_var for 1mu1p reconstructed variables initalized above

  // Initialize output trees
  // ... for nue intrinsic
  TTree* sel1e1p_nue_output_tree = new TTree("sel1e1p_nue_tree", "events selected as 1e1p BDT sideband from intrinsic nue overlay sample");
  sel1e1p_nue_output_tree->Branch("run",&run);
  sel1e1p_nue_output_tree->Branch("subrun",&subrun);
  sel1e1p_nue_output_tree->Branch("event",&event);
  sel1e1p_nue_output_tree->Branch("nu_energy_reco",&nu_energy_reco);
  sel1e1p_nue_output_tree->Branch("nu_pdg",&nu_pdg);
  sel1e1p_nue_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel1e1p_nue_output_tree->Branch("nu_L_true",&nu_L_true);
  sel1e1p_nue_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel1e1p_nue_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel1e1p_nue_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel1e1p_nue_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel1e1p_nue_output_tree->Branch("tune1_weight",&tune1_weight);
  sel1e1p_nue_output_tree->Branch("spline_weight",&spline_weight);
  sel1e1p_nue_output_tree->Branch("rootino_weight",&rootino_weight);
  sel1e1p_nue_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel1e1p_nue_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel1e1p_nue_output_tree->Branch("lee_weight",&lee_weight);
  sel1e1p_nue_output_tree->Branch("sys_weights",&sys_weights);
  // ... for nue bnb
  TTree* sel1e1p_bnb_output_tree = new TTree("sel1e1p_bnb_tree", "events selected as 1e1p BDT sideband from BNB overlay sample");
  sel1e1p_bnb_output_tree->Branch("run",&run);
  sel1e1p_bnb_output_tree->Branch("subrun",&subrun);
  sel1e1p_bnb_output_tree->Branch("event",&event);
  sel1e1p_bnb_output_tree->Branch("nu_energy_reco",&nu_energy_reco);
  sel1e1p_bnb_output_tree->Branch("nu_pdg",&nu_pdg);
  sel1e1p_bnb_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel1e1p_bnb_output_tree->Branch("nu_L_true",&nu_L_true);
  sel1e1p_bnb_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel1e1p_bnb_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel1e1p_bnb_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel1e1p_bnb_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel1e1p_bnb_output_tree->Branch("tune1_weight",&tune1_weight);
  sel1e1p_bnb_output_tree->Branch("spline_weight",&spline_weight);
  sel1e1p_bnb_output_tree->Branch("rootino_weight",&rootino_weight);
  sel1e1p_bnb_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel1e1p_bnb_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel1e1p_bnb_output_tree->Branch("lee_weight",&lee_weight);
  sel1e1p_bnb_output_tree->Branch("sys_weights",&sys_weights);
  // ... for nue dirt
  TTree* sel1e1p_dirt_output_tree = new TTree("sel1e1p_dirt_tree", "events selected as 1e1p BDT sideband from dirt overlay sample");
  sel1e1p_dirt_output_tree->Branch("run",&run);
  sel1e1p_dirt_output_tree->Branch("subrun",&subrun);
  sel1e1p_dirt_output_tree->Branch("event",&event);
  sel1e1p_dirt_output_tree->Branch("nu_energy_reco",&nu_energy_reco);
  sel1e1p_dirt_output_tree->Branch("nu_pdg",&nu_pdg);
  sel1e1p_dirt_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel1e1p_dirt_output_tree->Branch("nu_L_true",&nu_L_true);
  sel1e1p_dirt_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel1e1p_dirt_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel1e1p_dirt_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel1e1p_dirt_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel1e1p_dirt_output_tree->Branch("tune1_weight",&tune1_weight);
  sel1e1p_dirt_output_tree->Branch("spline_weight",&spline_weight);
  sel1e1p_dirt_output_tree->Branch("rootino_weight",&rootino_weight);
  sel1e1p_dirt_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel1e1p_dirt_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel1e1p_dirt_output_tree->Branch("lee_weight",&lee_weight);
  sel1e1p_dirt_output_tree->Branch("sys_weights",&sys_weights);
  // ... for nue extbnb
  TTree* sel1e1p_extbnb_output_tree = new TTree("sel1e1p_extbnb_tree", "events selected as 1e1p BDT sideband from EXTBNB data sample");
  sel1e1p_extbnb_output_tree->Branch("run",&run);
  sel1e1p_extbnb_output_tree->Branch("subrun",&subrun);
  sel1e1p_extbnb_output_tree->Branch("event",&event); 
  sel1e1p_extbnb_output_tree->Branch("tune1_weight",&tune1_weight);
  sel1e1p_extbnb_output_tree->Branch("nu_energy_reco",&nu_energy_reco);
  sel1e1p_extbnb_output_tree->Branch("spline_weight",&spline_weight);
  sel1e1p_extbnb_output_tree->Branch("rootino_weight",&rootino_weight);
  sel1e1p_extbnb_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel1e1p_extbnb_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel1e1p_extbnb_output_tree->Branch("lee_weight",&lee_weight);
  sel1e1p_extbnb_output_tree->Branch("sys_weights",&sys_weights);
  // ... for numu bnb
  TTree* sel1mu1p_bnb_output_tree = new TTree("sel1mu1p_bnb_tree", "events selected as 1mu1p from BNB overlay sample");
  sel1mu1p_bnb_output_tree->Branch("run",&run);
  sel1mu1p_bnb_output_tree->Branch("subrun",&subrun);
  sel1mu1p_bnb_output_tree->Branch("event",&event); 
  sel1mu1p_bnb_output_tree->Branch("nu_energy_reco",&reco_var[8]);
  sel1mu1p_bnb_output_tree->Branch("nu_pdg",&nu_pdg);
  sel1mu1p_bnb_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel1mu1p_bnb_output_tree->Branch("nu_L_true",&nu_L_true);
  sel1mu1p_bnb_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel1mu1p_bnb_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel1mu1p_bnb_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel1mu1p_bnb_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel1mu1p_bnb_output_tree->Branch("tune1_weight",&tune1_weight);
  sel1mu1p_bnb_output_tree->Branch("spline_weight",&spline_weight);
  sel1mu1p_bnb_output_tree->Branch("rootino_weight",&rootino_weight);
  sel1mu1p_bnb_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel1mu1p_bnb_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel1mu1p_bnb_output_tree->Branch("lee_weight",&lee_weight);
  sel1mu1p_bnb_output_tree->Branch("sys_weights",&sys_weights);
  sel1mu1p_bnb_output_tree->Branch("x_reco",&reco_var[0]);
  sel1mu1p_bnb_output_tree->Branch("y_reco",&reco_var[1]);
  sel1mu1p_bnb_output_tree->Branch("z_reco",&reco_var[2]);
  sel1mu1p_bnb_output_tree->Branch("eta_reco",&reco_var[3]);
  sel1mu1p_bnb_output_tree->Branch("openang_reco",&reco_var[4]);
  sel1mu1p_bnb_output_tree->Branch("sum_thetas_reco",&reco_var[5]);
  sel1mu1p_bnb_output_tree->Branch("sum_phis_reco",&reco_var[6]);
  sel1mu1p_bnb_output_tree->Branch("charge_near_trunk_reco",&reco_var[7]);
  sel1mu1p_bnb_output_tree->Branch("phiT_reco",&reco_var[9]);
  sel1mu1p_bnb_output_tree->Branch("alphaT_reco",&reco_var[10]);
  sel1mu1p_bnb_output_tree->Branch("pT_reco",&reco_var[11]);
  sel1mu1p_bnb_output_tree->Branch("pT_ratio_reco",&reco_var[12]);
  sel1mu1p_bnb_output_tree->Branch("Bjx_reco",&reco_var[13]);
  sel1mu1p_bnb_output_tree->Branch("Bjy_reco",&reco_var[14]);
  sel1mu1p_bnb_output_tree->Branch("Q2_reco",&reco_var[15]);
  sel1mu1p_bnb_output_tree->Branch("sph_reco",&reco_var[16]);
  sel1mu1p_bnb_output_tree->Branch("Q0_reco",&reco_var[17]);
  sel1mu1p_bnb_output_tree->Branch("Q3_reco",&reco_var[18]);
  sel1mu1p_bnb_output_tree->Branch("lepton_phi_reco",&reco_var[19]);
  sel1mu1p_bnb_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel1mu1p_bnb_output_tree->Branch("lepton_length_reco",&reco_var[21]);
  sel1mu1p_bnb_output_tree->Branch("lepton_KE_reco",&reco_var[22]);
  sel1mu1p_bnb_output_tree->Branch("lepton_cos_theta_reco",&reco_var[23]);
  sel1mu1p_bnb_output_tree->Branch("proton_phi_reco",&reco_var[24]);
  sel1mu1p_bnb_output_tree->Branch("proton_theta_reco",&reco_var[25]);
  sel1mu1p_bnb_output_tree->Branch("proton_length_reco",&reco_var[26]);
  sel1mu1p_bnb_output_tree->Branch("proton_KE_reco",&reco_var[27]);
  sel1mu1p_bnb_output_tree->Branch("proton_cos_theta_reco",&reco_var[28]);
  // ... for numu nue
  TTree* sel1mu1p_nue_output_tree = new TTree("sel1mu1p_nue_tree", "events selected as 1mu1p from intrinsic nue overlay sample");
  sel1mu1p_nue_output_tree->Branch("run",&run);
  sel1mu1p_nue_output_tree->Branch("subrun",&subrun);
  sel1mu1p_nue_output_tree->Branch("event",&event); 
  sel1mu1p_nue_output_tree->Branch("nu_energy_reco",&reco_var[8]);
  sel1mu1p_nue_output_tree->Branch("nu_pdg",&nu_pdg);
  sel1mu1p_nue_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel1mu1p_nue_output_tree->Branch("nu_L_true",&nu_L_true);
  sel1mu1p_nue_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel1mu1p_nue_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel1mu1p_nue_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel1mu1p_nue_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel1mu1p_nue_output_tree->Branch("tune1_weight",&tune1_weight);
  sel1mu1p_nue_output_tree->Branch("spline_weight",&spline_weight); 
  sel1mu1p_nue_output_tree->Branch("rootino_weight",&rootino_weight);
  sel1mu1p_nue_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel1mu1p_nue_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel1mu1p_nue_output_tree->Branch("lee_weight",&lee_weight);
  sel1mu1p_nue_output_tree->Branch("sys_weights",&sys_weights);
  sel1mu1p_nue_output_tree->Branch("x_reco",&reco_var[0]);
  sel1mu1p_nue_output_tree->Branch("y_reco",&reco_var[1]);
  sel1mu1p_nue_output_tree->Branch("z_reco",&reco_var[2]);
  sel1mu1p_nue_output_tree->Branch("eta_reco",&reco_var[3]);
  sel1mu1p_nue_output_tree->Branch("openang_reco",&reco_var[4]);
  sel1mu1p_nue_output_tree->Branch("sum_thetas_reco",&reco_var[5]);
  sel1mu1p_nue_output_tree->Branch("sum_phis_reco",&reco_var[6]);
  sel1mu1p_nue_output_tree->Branch("charge_near_trunk_reco",&reco_var[7]);
  sel1mu1p_nue_output_tree->Branch("phiT_reco",&reco_var[9]);
  sel1mu1p_nue_output_tree->Branch("alphaT_reco",&reco_var[10]);
  sel1mu1p_nue_output_tree->Branch("pT_reco",&reco_var[11]);
  sel1mu1p_nue_output_tree->Branch("pT_ratio_reco",&reco_var[12]);
  sel1mu1p_nue_output_tree->Branch("Bjx_reco",&reco_var[13]);
  sel1mu1p_nue_output_tree->Branch("Bjy_reco",&reco_var[14]);
  sel1mu1p_nue_output_tree->Branch("Q2_reco",&reco_var[15]);
  sel1mu1p_nue_output_tree->Branch("sph_reco",&reco_var[16]);
  sel1mu1p_nue_output_tree->Branch("Q0_reco",&reco_var[17]);
  sel1mu1p_nue_output_tree->Branch("Q3_reco",&reco_var[18]);
  sel1mu1p_nue_output_tree->Branch("lepton_phi_reco",&reco_var[19]);
  sel1mu1p_nue_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel1mu1p_nue_output_tree->Branch("lepton_length_reco",&reco_var[21]);
  sel1mu1p_nue_output_tree->Branch("lepton_KE_reco",&reco_var[22]);
  sel1mu1p_nue_output_tree->Branch("lepton_cos_theta_reco",&reco_var[23]);
  sel1mu1p_nue_output_tree->Branch("proton_phi_reco",&reco_var[24]);
  sel1mu1p_nue_output_tree->Branch("proton_theta_reco",&reco_var[25]);
  sel1mu1p_nue_output_tree->Branch("proton_length_reco",&reco_var[26]);
  sel1mu1p_nue_output_tree->Branch("proton_KE_reco",&reco_var[27]);
  sel1mu1p_nue_output_tree->Branch("proton_cos_theta_reco",&reco_var[28]);
  // ... for numu dirt
  TTree* sel1mu1p_dirt_output_tree = new TTree("sel1mu1p_dirt_tree", "events selected as 1mu1p from dirt overlay sample");
  sel1mu1p_dirt_output_tree->Branch("run",&run);
  sel1mu1p_dirt_output_tree->Branch("subrun",&subrun);
  sel1mu1p_dirt_output_tree->Branch("event",&event); 
  sel1mu1p_dirt_output_tree->Branch("nu_energy_reco",&reco_var[8]);
  sel1mu1p_dirt_output_tree->Branch("nu_pdg",&nu_pdg);
  sel1mu1p_dirt_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel1mu1p_dirt_output_tree->Branch("nu_L_true",&nu_L_true);
  sel1mu1p_dirt_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel1mu1p_dirt_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel1mu1p_dirt_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel1mu1p_dirt_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel1mu1p_dirt_output_tree->Branch("tune1_weight",&tune1_weight);
  sel1mu1p_dirt_output_tree->Branch("spline_weight",&spline_weight);
  sel1mu1p_dirt_output_tree->Branch("rootino_weight",&rootino_weight);
  sel1mu1p_dirt_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel1mu1p_dirt_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel1mu1p_dirt_output_tree->Branch("lee_weight",&lee_weight);
  sel1mu1p_dirt_output_tree->Branch("sys_weights",&sys_weights);
  sel1mu1p_dirt_output_tree->Branch("x_reco",&reco_var[0]);
  sel1mu1p_dirt_output_tree->Branch("y_reco",&reco_var[1]);
  sel1mu1p_dirt_output_tree->Branch("z_reco",&reco_var[2]);
  sel1mu1p_dirt_output_tree->Branch("eta_reco",&reco_var[3]);
  sel1mu1p_dirt_output_tree->Branch("openang_reco",&reco_var[4]);
  sel1mu1p_dirt_output_tree->Branch("sum_thetas_reco",&reco_var[5]);
  sel1mu1p_dirt_output_tree->Branch("sum_phis_reco",&reco_var[6]);
  sel1mu1p_dirt_output_tree->Branch("charge_near_trunk_reco",&reco_var[7]);
  sel1mu1p_dirt_output_tree->Branch("phiT_reco",&reco_var[9]);
  sel1mu1p_dirt_output_tree->Branch("alphaT_reco",&reco_var[10]);
  sel1mu1p_dirt_output_tree->Branch("pT_reco",&reco_var[11]);
  sel1mu1p_dirt_output_tree->Branch("pT_ratio_reco",&reco_var[12]);
  sel1mu1p_dirt_output_tree->Branch("Bjx_reco",&reco_var[13]);
  sel1mu1p_dirt_output_tree->Branch("Bjy_reco",&reco_var[14]);
  sel1mu1p_dirt_output_tree->Branch("Q2_reco",&reco_var[15]);
  sel1mu1p_dirt_output_tree->Branch("sph_reco",&reco_var[16]);
  sel1mu1p_dirt_output_tree->Branch("Q0_reco",&reco_var[17]);
  sel1mu1p_dirt_output_tree->Branch("Q3_reco",&reco_var[18]);
  sel1mu1p_dirt_output_tree->Branch("lepton_phi_reco",&reco_var[19]);
  sel1mu1p_dirt_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel1mu1p_dirt_output_tree->Branch("lepton_length_reco",&reco_var[21]);
  sel1mu1p_dirt_output_tree->Branch("lepton_KE_reco",&reco_var[22]);
  sel1mu1p_dirt_output_tree->Branch("lepton_cos_theta_reco",&reco_var[23]);
  sel1mu1p_dirt_output_tree->Branch("proton_phi_reco",&reco_var[24]);
  sel1mu1p_dirt_output_tree->Branch("proton_theta_reco",&reco_var[25]);
  sel1mu1p_dirt_output_tree->Branch("proton_length_reco",&reco_var[26]);
  sel1mu1p_dirt_output_tree->Branch("proton_KE_reco",&reco_var[27]);
  sel1mu1p_dirt_output_tree->Branch("proton_cos_theta_reco",&reco_var[28]);
  // ... for numu extbnb
  TTree* sel1mu1p_extbnb_output_tree = new TTree("sel1mu1p_extbnb_tree", "events selected as 1mu1p from EXTBNB data sample");
  sel1mu1p_extbnb_output_tree->Branch("run",&run);
  sel1mu1p_extbnb_output_tree->Branch("subrun",&subrun);
  sel1mu1p_extbnb_output_tree->Branch("event",&event); 
  sel1mu1p_extbnb_output_tree->Branch("nu_energy_reco",&reco_var[8]);
  sel1mu1p_extbnb_output_tree->Branch("tune1_weight",&tune1_weight);
  sel1mu1p_extbnb_output_tree->Branch("rootino_weight",&rootino_weight);
  sel1mu1p_extbnb_output_tree->Branch("spline_weight",&spline_weight);
  sel1mu1p_extbnb_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel1mu1p_extbnb_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel1mu1p_extbnb_output_tree->Branch("lee_weight",&lee_weight);
  sel1mu1p_extbnb_output_tree->Branch("sys_weights",&sys_weights);
  sel1mu1p_extbnb_output_tree->Branch("x_reco",&reco_var[0]);
  sel1mu1p_extbnb_output_tree->Branch("y_reco",&reco_var[1]);
  sel1mu1p_extbnb_output_tree->Branch("z_reco",&reco_var[2]);
  sel1mu1p_extbnb_output_tree->Branch("eta_reco",&reco_var[3]);
  sel1mu1p_extbnb_output_tree->Branch("openang_reco",&reco_var[4]);
  sel1mu1p_extbnb_output_tree->Branch("sum_thetas_reco",&reco_var[5]);
  sel1mu1p_extbnb_output_tree->Branch("sum_phis_reco",&reco_var[6]);
  sel1mu1p_extbnb_output_tree->Branch("charge_near_trunk_reco",&reco_var[7]);
  sel1mu1p_extbnb_output_tree->Branch("phiT_reco",&reco_var[9]);
  sel1mu1p_extbnb_output_tree->Branch("alphaT_reco",&reco_var[10]);
  sel1mu1p_extbnb_output_tree->Branch("pT_reco",&reco_var[11]);
  sel1mu1p_extbnb_output_tree->Branch("pT_ratio_reco",&reco_var[12]);
  sel1mu1p_extbnb_output_tree->Branch("Bjx_reco",&reco_var[13]);
  sel1mu1p_extbnb_output_tree->Branch("Bjy_reco",&reco_var[14]);
  sel1mu1p_extbnb_output_tree->Branch("Q2_reco",&reco_var[15]);
  sel1mu1p_extbnb_output_tree->Branch("sph_reco",&reco_var[16]);
  sel1mu1p_extbnb_output_tree->Branch("Q0_reco",&reco_var[17]);
  sel1mu1p_extbnb_output_tree->Branch("Q3_reco",&reco_var[18]);
  sel1mu1p_extbnb_output_tree->Branch("lepton_phi_reco",&reco_var[19]);
  sel1mu1p_extbnb_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel1mu1p_extbnb_output_tree->Branch("lepton_length_reco",&reco_var[21]);
  sel1mu1p_extbnb_output_tree->Branch("lepton_KE_reco",&reco_var[22]);
  sel1mu1p_extbnb_output_tree->Branch("lepton_cos_theta_reco",&reco_var[23]);
  sel1mu1p_extbnb_output_tree->Branch("proton_phi_reco",&reco_var[24]);
  sel1mu1p_extbnb_output_tree->Branch("proton_theta_reco",&reco_var[25]);
  sel1mu1p_extbnb_output_tree->Branch("proton_length_reco",&reco_var[26]);
  sel1mu1p_extbnb_output_tree->Branch("proton_KE_reco",&reco_var[27]);
  sel1mu1p_extbnb_output_tree->Branch("proton_cos_theta_reco",&reco_var[28]);
  // ... for numu bnb5e19 data
  TTree* sel1mu1p_5e19_output_tree = new TTree("sel1mu1p_5e19_tree", "events selected as 1mu1p from 5E19 data sample");
  sel1mu1p_5e19_output_tree->Branch("run",&run);
  sel1mu1p_5e19_output_tree->Branch("subrun",&subrun);
  sel1mu1p_5e19_output_tree->Branch("event",&event); 
  sel1mu1p_5e19_output_tree->Branch("nu_energy_reco",&nu_energy_reco);
  sel1mu1p_5e19_output_tree->Branch("tune1_weight",&tune1_weight);
  sel1mu1p_5e19_output_tree->Branch("spline_weight",&spline_weight);
  sel1mu1p_5e19_output_tree->Branch("rootino_weight",&rootino_weight);
  sel1mu1p_5e19_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel1mu1p_5e19_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel1mu1p_5e19_output_tree->Branch("lee_weight",&lee_weight);
  sel1mu1p_5e19_output_tree->Branch("sys_weights",&sys_weights);

  std::cout << "Initialized trees in output files" << std::endl;

  
  // * Get splines for Tune 1 reweighting * //
  
  // Open the v3.0.4 file, get the graphs
  TFile* v304_xsec_spline = new TFile(v304_xsec_spline_file, "READ");
  TGraph* xsec_mcc9_graph_numu = (TGraph*)v304_xsec_spline->Get("nu_mu_Ar40/qel_cc_n");
  TGraph* xsec_mcc9_graph_numubar = (TGraph*)v304_xsec_spline->Get("nu_mu_bar_Ar40/qel_cc_p");
  TGraph* xsec_mcc9_graph_nue = (TGraph*)v304_xsec_spline->Get("nu_e_Ar40/qel_cc_n");
  TGraph* xsec_mcc9_graph_nuebar = (TGraph*)v304_xsec_spline->Get("nu_e_bar_Ar40/qel_cc_p");
  // Open the Tune 1 file, get the graphs
  TFile* tune1_xsec_spline = new TFile(tune1_xsec_spline_file, "READ");
  TGraph* xsec_tune1_graph_numu = (TGraph*)tune1_xsec_spline->Get("nu_mu_Ar40/qel_cc_n");
  TGraph* xsec_tune1_graph_numubar = (TGraph*)tune1_xsec_spline->Get("nu_mu_bar_Ar40/qel_cc_p");
  TGraph* xsec_tune1_graph_nue = (TGraph*)tune1_xsec_spline->Get("nu_e_Ar40/qel_cc_n");
  TGraph* xsec_tune1_graph_nuebar = (TGraph*)tune1_xsec_spline->Get("nu_e_bar_Ar40/qel_cc_p");
  
  std::cout << "Got splines for Tune 1 reweighting" << std::endl;

  
  // * Fill trees in the output file * //
  
  // Read in nue arborist file
  TFile* nue_arborist = new TFile(nue_arborist_file, "READ");
  if ( !nue_arborist->IsOpen() ) std::cout << "Failed to open file " << nue_arborist_file << std::endl;
  
  // Get the desired tree
  TTree* nue_evweight_tree = (TTree*)nue_arborist->Get("arborist/eventweight_tree");
  // Get the desired branches
  nue_evweight_tree->SetBranchAddress("run",&run);
  nue_evweight_tree->SetBranchAddress("subrun",&subrun);
  nue_evweight_tree->SetBranchAddress("event",&event);
  nue_evweight_tree->SetBranchAddress("nu_pdg",&nu_pdg);
  nue_evweight_tree->SetBranchAddress("nu_energy_true",&nu_energy_true);
  nue_evweight_tree->SetBranchAddress("nu_interaction_ccnc",&nu_interaction_ccnc);
  nue_evweight_tree->SetBranchAddress("nu_interaction_mode",&nu_interaction_mode);
  nue_evweight_tree->SetBranchAddress("nu_interaction_type",&nu_interaction_type);
  nue_evweight_tree->SetBranchAddress("nu_target_pdg",&nu_target_pdg);
  nue_evweight_tree->SetBranchAddress("nu_L_true",&nu_L_true);
  nue_evweight_tree->SetBranchAddress("spline_weight",&spline_weight);
  nue_evweight_tree->SetBranchAddress("rootino_weight",&rootino_weight);
  nue_evweight_tree->SetBranchAddress("ub_tune_weight",&ub_tune_weight);
  nue_evweight_tree->SetBranchAddress("xsec_corr_weight",&xsec_corr_weight);
  nue_evweight_tree->SetBranchAddress("lee_weight",&lee_weight);
  nue_evweight_tree->SetBranchAddress("sys_weights",&sys_weights_ptr);
  
  // Loop over the nue tree
  for ( int i=0; i < nue_evweight_tree->GetEntries(); i++ ) {
    
    // Get the entry
    nue_evweight_tree->GetEntry(i);
    
    // Get the rse
    rse[0] = run;
    rse[1] = subrun;
    rse[2] = event;
    
    // If there's a corresponding entry in the 1e1p map...
    if ( sel1e1p_nue_event_map.find(rse) != sel1e1p_nue_event_map.end() ) {
      // Get the reconstructed energy
      nu_energy_reco = sel1e1p_nue_event_map.find(rse)->second;
      // Get the Tune 1 weight
      if ( nu_interaction_type == 1001 ) {
	if ( nu_pdg == 12 ) tune1_weight = xsec_tune1_graph_nue->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nue->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -12 ) tune1_weight = xsec_tune1_graph_nuebar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nuebar->Eval(nu_energy_true/1000.);
      }
      else tune1_weight = 1.;
      // If the ub_tune_weight is infinite or negative, set it to 1 -- per recommendation from Steven G.
      if ( std::isinf(ub_tune_weight) || ub_tune_weight < 0. ) {
	if ( std::isinf(ub_tune_weight) )
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1e1p_nue has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1e1p_nue has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the output tree
      sel1e1p_nue_output_tree->Fill();
    }

    // If there's a corresponding entry in the 1mu1p map...
    if ( sel1mu1p_nue_event_map.find(rse) != sel1mu1p_nue_event_map.end() ) {
      // Get the reconstructed variables
      reco_var = sel1mu1p_nue_event_map.find(rse)->second;
      // Get the Tune 1 weight
      if ( nu_interaction_type == 1001 ) {
	if ( nu_pdg == 12 ) tune1_weight = xsec_tune1_graph_nue->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nue->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -12 ) tune1_weight = xsec_tune1_graph_nuebar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nuebar->Eval(nu_energy_true/1000.);
      }
      else tune1_weight = 1.;
      // If the ub_tune_weight is infinite or negative, set it to 1 -- per recommendation from Steven G.
      if ( std::isinf(ub_tune_weight) || ub_tune_weight < 0. ) {
	if ( std::isinf(ub_tune_weight) )
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1mu1p_nue has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1mu1p_nue has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the output tree
      sel1mu1p_nue_output_tree->Fill();
    }
    
  }
  
  std::cout << "Filled nue trees" << std::endl;

  // Read in bnb arborist file
  TFile* bnb_arborist = new TFile(bnb_arborist_file, "READ");
  if ( !bnb_arborist->IsOpen() ) std::cout << "Failed to open file " << bnb_arborist_file << std::endl;
  
  // Get the desired tree
  TTree* bnb_evweight_tree = (TTree*)bnb_arborist->Get("arborist/eventweight_tree");
  bnb_evweight_tree->SetBranchAddress("run",&run);
  bnb_evweight_tree->SetBranchAddress("subrun",&subrun);
  bnb_evweight_tree->SetBranchAddress("event",&event);
  bnb_evweight_tree->SetBranchAddress("nu_pdg",&nu_pdg);
  bnb_evweight_tree->SetBranchAddress("nu_energy_true",&nu_energy_true);
  bnb_evweight_tree->SetBranchAddress("nu_interaction_ccnc",&nu_interaction_ccnc);
  bnb_evweight_tree->SetBranchAddress("nu_interaction_mode",&nu_interaction_mode);
  bnb_evweight_tree->SetBranchAddress("nu_interaction_type",&nu_interaction_type);
  bnb_evweight_tree->SetBranchAddress("nu_target_pdg",&nu_target_pdg);
  bnb_evweight_tree->SetBranchAddress("nu_L_true",&nu_L_true);
  bnb_evweight_tree->SetBranchAddress("spline_weight",&spline_weight);
  bnb_evweight_tree->SetBranchAddress("rootino_weight",&rootino_weight);
  bnb_evweight_tree->SetBranchAddress("ub_tune_weight",&ub_tune_weight);
  bnb_evweight_tree->SetBranchAddress("xsec_corr_weight",&xsec_corr_weight);
  bnb_evweight_tree->SetBranchAddress("lee_weight",&lee_weight);
  bnb_evweight_tree->SetBranchAddress("sys_weights",&sys_weights_ptr);
  
  // Loop over the bnb tree
  for ( int i=0; i < bnb_evweight_tree->GetEntries(); i++ ) {
    
    // Get the entry
    bnb_evweight_tree->GetEntry(i);
    
    // Get the rse
    rse[0] = run;
    rse[1] = subrun;
    rse[2] = event;
    
    // If there's a corresponding entry in the 1e1p map...
    if ( sel1e1p_bnb_event_map.find(rse) != sel1e1p_bnb_event_map.end() ) {
      // Get the reconstructed energy
      nu_energy_reco = sel1e1p_bnb_event_map.find(rse)->second;
      // Get the Tune 1 weight
      if ( nu_interaction_type == 1001 ) {
	if ( nu_pdg == 14 ) tune1_weight = xsec_tune1_graph_numu->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_numu->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -14 ) tune1_weight = xsec_tune1_graph_numubar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_numubar->Eval(nu_energy_true/1000.);
	if ( nu_pdg == 12 ) tune1_weight = xsec_tune1_graph_nue->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nue->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -12 ) tune1_weight = xsec_tune1_graph_nuebar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nuebar->Eval(nu_energy_true/1000.);
      }
      else tune1_weight = 1.;
      // If the ub_tune_weight is infinite or negative, set it to 1 -- per recommendation from Steven G.
      if ( std::isinf(ub_tune_weight) || ub_tune_weight < 0. ) {
	if ( std::isinf(ub_tune_weight) )
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1e1p_bnb has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1e1p_bnb has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the tree
      sel1e1p_bnb_output_tree->Fill();
    }
    
    // If there's a corresponding entry in the 1mu1p map...
    if ( sel1mu1p_bnb_event_map.find(rse) != sel1mu1p_bnb_event_map.end() ) {
      // Get the reconstructed variables
      reco_var = sel1mu1p_bnb_event_map.find(rse)->second;
      // Get the Tune 1 weight
      if ( nu_interaction_type == 1001 ) {
	if ( nu_pdg == 14 ) tune1_weight = xsec_tune1_graph_numu->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_numu->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -14 ) tune1_weight = xsec_tune1_graph_numubar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_numubar->Eval(nu_energy_true/1000.);
	if ( nu_pdg == 12 ) tune1_weight = xsec_tune1_graph_nue->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nue->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -12 ) tune1_weight = xsec_tune1_graph_nuebar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nuebar->Eval(nu_energy_true/1000.);
      }
      else tune1_weight = 1.;
      // If the ub_tune_weight is infinite or negative, set it to 1 -- per recommendation from Steven G.
      if ( std::isinf(ub_tune_weight) || ub_tune_weight < 0. ) {
	if ( std::isinf(ub_tune_weight) )
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1mu1p_bnb has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1mu1p_bnb has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the tree
      sel1mu1p_bnb_output_tree->Fill();
    }
    
  }
  
  std::cout << "Filled bnb mc trees" << std::endl;

  // Read in dirt arborist file
  TFile* dirt_arborist = new TFile(dirt_arborist_file, "READ");
  if ( !dirt_arborist->IsOpen() ) std::cout << "Failed to open file " << dirt_arborist_file << std::endl;
  
  // Get the desired tree
  TTree* dirt_evweight_tree = (TTree*)dirt_arborist->Get("arborist/eventweight_tree");
  dirt_evweight_tree->SetBranchAddress("run",&run);
  dirt_evweight_tree->SetBranchAddress("subrun",&subrun);
  dirt_evweight_tree->SetBranchAddress("event",&event);
  dirt_evweight_tree->SetBranchAddress("nu_pdg",&nu_pdg);
  dirt_evweight_tree->SetBranchAddress("nu_energy_true",&nu_energy_true);
  dirt_evweight_tree->SetBranchAddress("nu_interaction_ccnc",&nu_interaction_ccnc);
  dirt_evweight_tree->SetBranchAddress("nu_interaction_mode",&nu_interaction_mode);
  dirt_evweight_tree->SetBranchAddress("nu_interaction_type",&nu_interaction_type);
  dirt_evweight_tree->SetBranchAddress("nu_target_pdg",&nu_target_pdg);
  dirt_evweight_tree->SetBranchAddress("nu_L_true",&nu_L_true);
  dirt_evweight_tree->SetBranchAddress("spline_weight",&spline_weight);
  dirt_evweight_tree->SetBranchAddress("rootino_weight",&rootino_weight);
  dirt_evweight_tree->SetBranchAddress("ub_tune_weight",&ub_tune_weight);
  dirt_evweight_tree->SetBranchAddress("xsec_corr_weight",&xsec_corr_weight);
  dirt_evweight_tree->SetBranchAddress("lee_weight",&lee_weight);
  dirt_evweight_tree->SetBranchAddress("sys_weights",&sys_weights_ptr);
  
  // Loop over the dirt tree
  for ( int i=0; i < dirt_evweight_tree->GetEntries(); i++ ) {
    
    // Get the entry
    dirt_evweight_tree->GetEntry(i);
    
    // Get the rse
    rse[0] = run;
    rse[1] = subrun;
    rse[2] = event;
    
    // If there's a corresponding entry in the 1e1p map...
    if ( sel1e1p_dirt_event_map.find(rse) != sel1e1p_dirt_event_map.end() ) {
      // Get the reconstructed energy
      nu_energy_reco = sel1e1p_dirt_event_map.find(rse)->second;
      // Get the Tune 1 weight
      if ( nu_interaction_type == 1001 ) {
	if ( nu_pdg == 14 ) tune1_weight = xsec_tune1_graph_numu->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_numu->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -14 ) tune1_weight = xsec_tune1_graph_numubar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_numubar->Eval(nu_energy_true/1000.);
	if ( nu_pdg == 12 ) tune1_weight = xsec_tune1_graph_nue->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nue->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -12 ) tune1_weight = xsec_tune1_graph_nuebar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nuebar->Eval(nu_energy_true/1000.);
      }
      else tune1_weight = 1.;
      // If the ub_tune_weight is infinite or negative, set it to 1 -- per recommendation from Steven G.
      if ( std::isinf(ub_tune_weight) || ub_tune_weight < 0. ) {
	if ( std::isinf(ub_tune_weight) )
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1e1p_dirt has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1e1p_dirt has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the tree
      sel1e1p_dirt_output_tree->Fill();
    }
    
    // If there's a corresponding entry in the 1mu1p map...
    if ( sel1mu1p_dirt_event_map.find(rse) != sel1mu1p_dirt_event_map.end() ) {
      // Get the reconstructed variables
      reco_var = sel1mu1p_dirt_event_map.find(rse)->second;
      // Get the Tune 1 weight
      if ( nu_interaction_type == 1001 ) {
	if ( nu_pdg == 14 ) tune1_weight = xsec_tune1_graph_numu->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_numu->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -14 ) tune1_weight = xsec_tune1_graph_numubar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_numubar->Eval(nu_energy_true/1000.);
	if ( nu_pdg == 12 ) tune1_weight = xsec_tune1_graph_nue->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nue->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -12 ) tune1_weight = xsec_tune1_graph_nuebar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nuebar->Eval(nu_energy_true/1000.);
      }
      else tune1_weight = 1.;
      // If the ub_tune_weight is infinite or negative, set it to 1 -- per recommendation from Steven G.
      if ( std::isinf(ub_tune_weight) || ub_tune_weight < 0. ) {
	if ( std::isinf(ub_tune_weight) )
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1mu1p_dirt has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel1mu1p_dirt has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the tree
      sel1mu1p_dirt_output_tree->Fill();
    }
    
  }
  
  std::cout << "Filled dirt trees" << std::endl;
  
  
  // For data samples, set dummy values for various weights
  spline_weight = 1.;
  rootino_weight = 1.;
  ub_tune_weight = 1.;
  xsec_corr_weight = 1.;
  tune1_weight = 1.;
  lee_weight = 0.;
  for ( auto weight_it = sys_weights.begin(); weight_it != sys_weights.end(); weight_it++ ) {
    for ( unsigned int i=0; i < weight_it->second.size(); i++ ) {
      weight_it->second[i] = 1.;
    }
  }
  
  // Fill 1e1p extbnb tree
  for ( auto event_it = sel1e1p_extbnb_event_map.begin(); event_it != sel1e1p_extbnb_event_map.end(); event_it++ ) {
    // Get the variables to be put into the tree                                                                                                                                 
    run = event_it->first[0];
    subrun = event_it->first[1];
    event = event_it->first[2];
    nu_energy_reco = event_it->second;
    // Fill the tree                                                                                                                                                             
    sel1e1p_extbnb_output_tree->Fill();
  }
  
  // Fill 1mu1p extbnb tree
  for ( auto event_it = sel1mu1p_extbnb_event_map.begin(); event_it != sel1mu1p_extbnb_event_map.end(); event_it++ ) {
    // Get the variables to be put into the tree                                                                                                                                 
    run = event_it->first[0];
    subrun = event_it->first[1];
    event = event_it->first[2];
    reco_var = event_it->second;
    // Fill the tree                                                                                                                                                             
    sel1mu1p_extbnb_output_tree->Fill();
  }
  
  // Fill 1mu1p bnb5e19 tree
  for ( auto event_it = sel1mu1p_5e19_event_map.begin(); event_it != sel1mu1p_5e19_event_map.end(); event_it++ ) {
    // Get the variables to be put into the tree                                                                                                                                 
    run = event_it->first[0];
    subrun = event_it->first[1];
    event = event_it->first[2];
    nu_energy_reco = event_it->second;
    // Fill the tree                                                                                                                                                             
    sel1mu1p_5e19_output_tree->Fill();
  }

  std::cout << "Filled data trees" << std::endl;


  // * Conclusion * //
  
  // Count events in all output trees and make sure we did this right
  // ... for 1e1p
  int n_sel1e1p_nue_in     = sel1e1p_nue_event_map.size();
  int n_sel1e1p_nue_out    = sel1e1p_nue_output_tree->GetEntries();
  int n_sel1e1p_bnb_in     = sel1e1p_bnb_event_map.size();
  int n_sel1e1p_bnb_out    = sel1e1p_bnb_output_tree->GetEntries();
  int n_sel1e1p_dirt_in    = sel1e1p_dirt_event_map.size();
  int n_sel1e1p_dirt_out   = sel1e1p_dirt_output_tree->GetEntries();
  int n_sel1e1p_extbnb_in  = sel1e1p_extbnb_event_map.size();
  int n_sel1e1p_extbnb_out = sel1e1p_extbnb_output_tree->GetEntries();
  // ... for 1mu1p
  int n_sel1mu1p_bnb_in     = sel1mu1p_bnb_event_map.size();
  int n_sel1mu1p_bnb_out    = sel1mu1p_bnb_output_tree->GetEntries();
  int n_sel1mu1p_nue_in     = sel1mu1p_nue_event_map.size();
  int n_sel1mu1p_nue_out    = sel1mu1p_nue_output_tree->GetEntries();
  int n_sel1mu1p_dirt_in    = sel1mu1p_dirt_event_map.size();
  int n_sel1mu1p_dirt_out   = sel1mu1p_dirt_output_tree->GetEntries();
  int n_sel1mu1p_extbnb_in  = sel1mu1p_extbnb_event_map.size();
  int n_sel1mu1p_extbnb_out = sel1mu1p_extbnb_output_tree->GetEntries();
  // ... for 1mu1p bnb5e19
  int n_sel1mu1p_5e19_in  = sel1mu1p_5e19_event_map.size();
  int n_sel1mu1p_5e19_out = sel1mu1p_5e19_output_tree->GetEntries();
  
  // Print information
  // ... for 1e1p
  std::cout << "Selected 1e1p event count information" << std::endl;
  std::cout << "  intrinsic nue overlay events: " << n_sel1e1p_nue_in << " events in, " << n_sel1e1p_nue_out << " events out" << std::endl;
  std::cout << "  BNB overlay events: " << n_sel1e1p_bnb_in << " events in, " << n_sel1e1p_bnb_out << " events out" << std::endl;
  std::cout << "  dirt overlay events: " << n_sel1e1p_dirt_in << " events in, " << n_sel1e1p_dirt_out << " events out" << std::endl;
  std::cout << "  EXTBNB data events: " << n_sel1e1p_extbnb_in << " events in, " << n_sel1e1p_extbnb_out << " events out" << std::endl;
  if ( (n_sel1e1p_nue_in == n_sel1e1p_nue_out) and (n_sel1e1p_bnb_in == n_sel1e1p_bnb_out) and
       (n_sel1e1p_dirt_in == n_sel1e1p_dirt_out) and (n_sel1e1p_extbnb_in == n_sel1e1p_extbnb_out) ) {
    std::cout << "These values are consistent!" << std::endl;
  }
  else { std::cout << "Warning: Number of events in 1e1p output trees is inconsistent with input lists" << std::endl; }
  // ... for 1mu1p
  std::cout << "Selected 1mu1p event count information" << std::endl;
  std::cout << "  BNB overlay events: " << n_sel1mu1p_bnb_in << " events in, " << n_sel1mu1p_bnb_out << " events out" << std::endl;
  std::cout << "  intrinsic nue overlay events: " << n_sel1mu1p_nue_in << " events in, " << n_sel1mu1p_nue_out << " events out" << std::endl;
  std::cout << "  dirt overlay events: " << n_sel1mu1p_dirt_in << " events in, " << n_sel1mu1p_dirt_out << " events out" << std::endl;
  std::cout << "  EXTBNB data events: " << n_sel1mu1p_extbnb_in << " events in, " << n_sel1mu1p_extbnb_out << " events out" << std::endl;
  if ( (n_sel1mu1p_bnb_in == n_sel1mu1p_bnb_out) and (n_sel1mu1p_dirt_in == n_sel1mu1p_dirt_out) and (n_sel1mu1p_extbnb_in == n_sel1mu1p_extbnb_out) ) {
    if ( n_sel1mu1p_nue_in == n_sel1mu1p_nue_out ) std::cout << "These values are consistent!" << std::endl;
  }
  else { std::cout << "Warning: Number of events in 1mu1p output trees is inconsistent with input list" << std::endl; }
  // ... for 1mu1p bnb5e19
  std::cout << "Selected 1mu1p in BNB 5e19 event count information" << std::endl;
  std::cout << "  BNB 5e19 events: " << n_sel1mu1p_5e19_in << " events in, " << n_sel1mu1p_5e19_out << " events out" << std::endl;
  if ( n_sel1mu1p_5e19_in == n_sel1mu1p_5e19_out ) { std::cout << "These values are consistent!" << std::endl; }
  else { std::cout << "Warning: Number of events in BNB 5e19  output tree is inconsistent with input list" << std::endl; }
  
  // Write and close the output file
  out->Write();
  out->Close();
  
  // El Fin
  std::cout << "Wrote output file " << output_file << std::endl;
  std::cout << "Done!" << std::endl;
  return 0;
  
}
