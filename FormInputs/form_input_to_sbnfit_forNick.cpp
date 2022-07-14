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
  const char* sel_bnb_eventlist_file      = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/sel1e1p_FinalSelection_bnb_run3_eventlist.txt";
  const char* sel_bnb_lowE_eventlist_file = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/empty.txt";
  const char* sel_nue_eventlist_file      = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/sel1e1p_FinalSelection_nue_run3_eventlist.txt";
  const char* sel_nue_lowE_eventlist_file = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/empty.txt";
  const char* sel_ncpi0_eventlist_file    = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/sel1e1p_FinalSelection_ncpi0_run3_eventlist.txt";
  const char* sel_ccpi0_eventlist_file    = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/sel1e1p_FinalSelection_ccpi0_run3_eventlist.txt";
  const char* sel_dirt_eventlist_file     = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/empty.txt";
  const char* sel_extbnb_eventlist_file   = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/sel1e1p_FinalSelection_extbnb_run3_eventlist.txt";
  const char* sel_fullosc_eventlist_file   = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/sel1e1p_FinalSelection_fullosc_run3_eventlist.txt";
  // ... for data event list
  const char* sel_data_eventlist_file = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/empty.txt";
  // ... for Tune 1 weightsi
  const char* v304_xsec_spline_file  = "/uboone/data/users/yatesla/othersys_mcc9/for_Tune1_weights/xsec_graphs_mcc9_v304.root";
  const char* v304a_xsec_spline_file = "/uboone/data/users/yatesla/othersys_mcc9/for_Tune1_weights/xsec_graphs_mcc9_v304a.root";
  const char* tune1_xsec_spline_file = "/uboone/data/users/yatesla/othersys_mcc9/for_Tune1_weights/xsec_graphs_tune1.root";
  // ... for arborist files
  const char* bnb_arborist_file       = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v55_Jun20_withExtraGENIE_bnb_nu_run3.root";
  const char* bnb_lowE_arborist_file  = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v55_Jun20_withExtraGENIE_bnb_nu_lowE_run3.root";
  const char* nue_arborist_file       = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v55_Jun20_withExtraGENIE_intrinsic_nue_run3.root";
  const char* nue_lowE_arborist_file  = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v55_Jun20_withExtraGENIE_intrinsic_nue_lowE_run3.root";
  const char* ncpi0_arborist_file     = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v55_Jun20_withExtraGENIE_nc_pi0_run3.root";
  const char* ccpi0_arborist_file     = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v55_Jun20_withExtraGENIE_cc_pi0_run3.root";
  const char* dirt_arborist_file      = "/uboone/data/users/yatesla/othersys_mcc9/arborist/arborist_v55_Jun20_withExtraGENIE_dirt_nu_run3.root";
  // ... for output
  const char* output_file = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit/input_to_sbnfit_v55_Oct05_withExtraGENIE_1e1p_FinalSelection_run3.root";

  // Initalize the number of reco variables
  const int N_var = 93;  // 35 reco variables + true energy + ccnc flag + interaction type + scedr + "POT weight" + CV weight + LEE weight + pi0 flag (x2) + Run number
                         // + cut flag + BDT score information (46 variables) + pi0 weight
  //const int N_var = 2;

  
  //** Read in lists of selected events and their reconstructed energies **//
  
  // Initialize maps to hold selected events
  std::map<std::array<int,3>,std::array<double,N_var>> sel_bnb_event_map;
  std::map<std::array<int,3>,std::array<double,N_var>> sel_bnb_lowE_event_map;
  std::map<std::array<int,3>,std::array<double,N_var>> sel_nue_event_map;
  std::map<std::array<int,3>,std::array<double,N_var>> sel_nue_lowE_event_map;
  std::map<std::array<int,3>,std::array<double,N_var>> sel_ncpi0_event_map;
  std::map<std::array<int,3>,std::array<double,N_var>> sel_ccpi0_event_map;
  std::map<std::array<int,3>,std::array<double,N_var>> sel_dirt_event_map;
  std::map<std::array<int,3>,std::array<double,N_var>> sel_extbnb_event_map;
  std::map<std::array<int,3>,double> sel_data_event_map;
  
  // Initialize variables to hold intermediate values
  std::string line;
  std::vector<std::string> split_line;
  std::array<int,3> rse;
  double reco_energy;
  std::array<double,N_var> reco_var;

  // Read in the bnb nu event list
  std::ifstream sel_bnb_eventlist(sel_bnb_eventlist_file);
  if ( sel_bnb_eventlist.is_open() ) {
    while ( getline(sel_bnb_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      //std::cout << rse[0] << "," << rse[1] << "," << rse[2] << std::endl;
      for ( int i=0; i<N_var; i++ ) {
	//std::cout << i << ": " << split_line[i+4] << std::endl;
	reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
      }
      sel_bnb_event_map.emplace(rse,reco_var);
    }
    sel_bnb_eventlist.close();
  }
  std::cout << "Read in " << sel_bnb_event_map.size() << " selected events from the BNB nu overlay sample" << std::endl;

  // Read in the bnb nu lowE event list
  std::ifstream sel_bnb_lowE_eventlist(sel_bnb_lowE_eventlist_file);
  if ( sel_bnb_lowE_eventlist.is_open() ) {
    while ( getline(sel_bnb_lowE_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<N_var; i++ ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
      sel_bnb_lowE_event_map.emplace(rse,reco_var);
    }
    sel_bnb_lowE_eventlist.close();
  }
  std::cout << "Read in " << sel_bnb_lowE_event_map.size() << " selected events from the BNB nu lowE overlay sample" << std::endl;

  // Read in the nue intrinsic event list
  std::ifstream sel_nue_eventlist(sel_nue_eventlist_file);
  if ( sel_nue_eventlist.is_open() ) {
    while ( getline(sel_nue_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      //std::cout << rse[0] << "," << rse[1] << "," << rse[2] << std::endl;
      for ( int i=0; i<N_var; i++ ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
      sel_nue_event_map.emplace(rse,reco_var);
    }
    sel_nue_eventlist.close();
  }
  std::cout << "Read in " << sel_nue_event_map.size() << " selected events from the intinsic nue overlay sample" << std::endl;

  // Read in the nue intrinsic lowE event list
  std::ifstream sel_nue_lowE_eventlist(sel_nue_lowE_eventlist_file);
  if ( sel_nue_lowE_eventlist.is_open() ) {
    while ( getline(sel_nue_lowE_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<N_var; i++ ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
      sel_nue_lowE_event_map.emplace(rse,reco_var);
    }
    sel_nue_lowE_eventlist.close();
  }
  std::cout << "Read in " << sel_nue_lowE_event_map.size() << " selected events from the intinsic nue lowE overlay sample" << std::endl;

  // Read in the NC pi0 event list
  std::ifstream sel_ncpi0_eventlist(sel_ncpi0_eventlist_file);
  if ( sel_ncpi0_eventlist.is_open() ) {
    while ( getline(sel_ncpi0_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<N_var; i++ ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
      sel_ncpi0_event_map.emplace(rse,reco_var);
    }
    sel_ncpi0_eventlist.close();
  }
  std::cout << "Read in " << sel_ncpi0_event_map.size() << " selected events from the NC pi0 overlay sample" << std::endl;

  // Read in the CC pi0 event list
  std::ifstream sel_ccpi0_eventlist(sel_ccpi0_eventlist_file);
  if ( sel_ccpi0_eventlist.is_open() ) {
    while ( getline(sel_ccpi0_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      for ( int i=0; i<N_var; i++ ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
      sel_ccpi0_event_map.emplace(rse,reco_var);
    }
    sel_ccpi0_eventlist.close();
  }
  std::cout << "Read in " << sel_ccpi0_event_map.size() << " selected events from the CC pi0 overlay sample" << std::endl;
    
  // Read in the dirt event list
  std::ifstream sel_dirt_eventlist(sel_dirt_eventlist_file);
  reco_var[N_var-1] = 1.0;  // set pi0 weights to 1.0
  if ( sel_dirt_eventlist.is_open() ) {
    while ( getline(sel_dirt_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      // -1 for no pi0 weights
      for ( int i=0; i<N_var-1; i++ ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
      sel_dirt_event_map.emplace(rse,reco_var);
    }
    sel_dirt_eventlist.close();
  }
  std::cout << "Read in " << sel_dirt_event_map.size() << " selected events from the dirt overlay sample" << std::endl;
  
  // Read in the extbnb event list
  std::ifstream sel_extbnb_eventlist(sel_extbnb_eventlist_file);
  reco_var[N_var-1] = 1.0;  // set pi0 weights to 1.0
  if ( sel_extbnb_eventlist.is_open() ) {
    while ( getline(sel_extbnb_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      // -1 for no pi0 weights
      for ( int i=0; i<N_var-1; i++ ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
      sel_extbnb_event_map.emplace(rse,reco_var);
    }
    sel_extbnb_eventlist.close();
  }
  std::cout << "Read in " << sel_extbnb_event_map.size() << " selected events from the EXTBNB data sample" << std::endl;
  
  // Read in the bnb data event list
  std::ifstream sel_data_eventlist(sel_data_eventlist_file);
  if ( sel_data_eventlist.is_open() ) {
    while ( getline(sel_data_eventlist, line) ) {
      boost::algorithm::split(split_line, line, boost::is_any_of(","));
      rse[0] = std::stoi(split_line[0]);
      rse[1] = std::stoi(split_line[1]);
      rse[2] = std::stoi(split_line[2]);
      reco_energy = std::stof(split_line[3]);
      sel_data_event_map.emplace(rse,reco_energy);
    }
    sel_data_eventlist.close();
  }
  std::cout << "Read in " << sel_data_event_map.size() << " selected events from the BNB data sample" << std::endl;
  
  
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

  // Initialize output trees
  // ... for bnb nu
  TTree* sel_bnb_output_tree = new TTree("sel_bnb_tree", "events selected from BNB nu overlay sample");
  sel_bnb_output_tree->Branch("run",&run);
  sel_bnb_output_tree->Branch("subrun",&subrun);
  sel_bnb_output_tree->Branch("event",&event);
  sel_bnb_output_tree->Branch("nu_pdg",&nu_pdg);
  sel_bnb_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel_bnb_output_tree->Branch("nu_L_true",&nu_L_true);
  sel_bnb_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel_bnb_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel_bnb_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel_bnb_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel_bnb_output_tree->Branch("tune1_weight",&tune1_weight);
  sel_bnb_output_tree->Branch("spline_weight",&spline_weight);
  sel_bnb_output_tree->Branch("rootino_weight",&rootino_weight);
  sel_bnb_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel_bnb_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel_bnb_output_tree->Branch("lee_weight",&lee_weight);
  sel_bnb_output_tree->Branch("sys_weights",&sys_weights);
  sel_bnb_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_bnb_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_bnb_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_bnb_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_bnb_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_bnb_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_bnb_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_bnb_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_bnb_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_bnb_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_bnb_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_bnb_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_bnb_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_bnb_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_bnb_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_bnb_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_bnb_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_bnb_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_bnb_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_bnb_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_bnb_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_bnb_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_bnb_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_bnb_output_tree->Branch("x_reco",&reco_var[23]);
  sel_bnb_output_tree->Branch("y_reco",&reco_var[24]);
  sel_bnb_output_tree->Branch("z_reco",&reco_var[25]);
  sel_bnb_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_bnb_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_bnb_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_bnb_output_tree->Branch("shr_charge_ratio_reco",&reco_var[29]);
  sel_bnb_output_tree->Branch("shr_consistency_reco",&reco_var[30]);
  sel_bnb_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[31]);
  sel_bnb_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[32]);
  sel_bnb_output_tree->Branch("proton_length_reco",&reco_var[33]);
  sel_bnb_output_tree->Branch("lepton_length_reco",&reco_var[34]);
  // skip reco_var[35], which is Nick's copy of the true neutrino energy
  // skip reco_var[36], which is Nick's copy of the ccnc flag
  // skip reco_var[37], which is Nick's copy of the interaction type
  sel_bnb_output_tree->Branch("dllee_scedr",&reco_var[38]);
  sel_bnb_output_tree->Branch("dllee_pot_weight",&reco_var[39]);
  // skip reco_var[40], which is Nick's copy of the CV tune weight
  // skip reco_var[41], which is Nick's copy of the LEE weight
  sel_bnb_output_tree->Branch("dllee_pi0_flag_new",&reco_var[42]);
  sel_bnb_output_tree->Branch("dllee_pi0_flag_old",&reco_var[43]);
  // skip reco_var[44], which is the Run number
  sel_bnb_output_tree->Branch("dllee_cut_level",&reco_var[45]);
  sel_bnb_output_tree->Branch("dllee_nBDTs",&reco_var[46]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[47]);
  sel_bnb_output_tree->Branch("dllee_bdt_score00",&reco_var[48]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[49]);
  sel_bnb_output_tree->Branch("dllee_bdt_score01",&reco_var[50]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[51]);
  sel_bnb_output_tree->Branch("dllee_bdt_score02",&reco_var[52]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[53]);
  sel_bnb_output_tree->Branch("dllee_bdt_score03",&reco_var[54]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[55]);
  sel_bnb_output_tree->Branch("dllee_bdt_score04",&reco_var[56]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[57]);
  sel_bnb_output_tree->Branch("dllee_bdt_score05",&reco_var[58]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[59]);
  sel_bnb_output_tree->Branch("dllee_bdt_score06",&reco_var[60]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[61]);
  sel_bnb_output_tree->Branch("dllee_bdt_score07",&reco_var[62]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[63]);
  sel_bnb_output_tree->Branch("dllee_bdt_score08",&reco_var[64]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[65]);
  sel_bnb_output_tree->Branch("dllee_bdt_score09",&reco_var[66]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight10",&reco_var[67]);
  sel_bnb_output_tree->Branch("dllee_bdt_score10",&reco_var[68]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight11",&reco_var[69]);
  sel_bnb_output_tree->Branch("dllee_bdt_score11",&reco_var[70]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight12",&reco_var[71]);
  sel_bnb_output_tree->Branch("dllee_bdt_score12",&reco_var[72]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight13",&reco_var[73]);
  sel_bnb_output_tree->Branch("dllee_bdt_score13",&reco_var[74]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight14",&reco_var[75]);
  sel_bnb_output_tree->Branch("dllee_bdt_score14",&reco_var[76]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight15",&reco_var[77]);
  sel_bnb_output_tree->Branch("dllee_bdt_score15",&reco_var[78]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight16",&reco_var[79]);
  sel_bnb_output_tree->Branch("dllee_bdt_score16",&reco_var[80]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight17",&reco_var[81]);
  sel_bnb_output_tree->Branch("dllee_bdt_score17",&reco_var[82]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight18",&reco_var[83]);
  sel_bnb_output_tree->Branch("dllee_bdt_score18",&reco_var[84]);
  sel_bnb_output_tree->Branch("dllee_bdt_tvweight19",&reco_var[85]);
  sel_bnb_output_tree->Branch("dllee_bdt_score19",&reco_var[86]);
  sel_bnb_output_tree->Branch("dllee_bdt_score_avg",&reco_var[87]);
  sel_bnb_output_tree->Branch("dllee_bdt_score_median",&reco_var[88]);
  sel_bnb_output_tree->Branch("dllee_bdt_score_max",&reco_var[89]);
  // skip reco_var[90] and reco_var[91], which are not clear to me
  sel_bnb_output_tree->Branch("dllee_pi0_weight",&reco_var[92]);
  // ... for bnb nu lowE
  TTree* sel_bnb_lowE_output_tree = new TTree("sel_bnb_lowE_tree", "events selected from BNB nu lowE overlay sample");
  sel_bnb_lowE_output_tree->Branch("run",&run);
  sel_bnb_lowE_output_tree->Branch("subrun",&subrun);
  sel_bnb_lowE_output_tree->Branch("event",&event);
  sel_bnb_lowE_output_tree->Branch("nu_pdg",&nu_pdg);
  sel_bnb_lowE_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel_bnb_lowE_output_tree->Branch("nu_L_true",&nu_L_true);
  sel_bnb_lowE_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel_bnb_lowE_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel_bnb_lowE_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel_bnb_lowE_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel_bnb_lowE_output_tree->Branch("tune1_weight",&tune1_weight);
  sel_bnb_lowE_output_tree->Branch("spline_weight",&spline_weight);
  sel_bnb_lowE_output_tree->Branch("rootino_weight",&rootino_weight);
  sel_bnb_lowE_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel_bnb_lowE_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel_bnb_lowE_output_tree->Branch("lee_weight",&lee_weight);
  sel_bnb_lowE_output_tree->Branch("sys_weights",&sys_weights);
  sel_bnb_lowE_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_bnb_lowE_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_bnb_lowE_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_bnb_lowE_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_bnb_lowE_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_bnb_lowE_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_bnb_lowE_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_bnb_lowE_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_bnb_lowE_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_bnb_lowE_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_bnb_lowE_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_bnb_lowE_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_bnb_lowE_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_bnb_lowE_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_bnb_lowE_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_bnb_lowE_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_bnb_lowE_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_bnb_lowE_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_bnb_lowE_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_bnb_lowE_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_bnb_lowE_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_bnb_lowE_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_bnb_lowE_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_bnb_lowE_output_tree->Branch("x_reco",&reco_var[23]);
  sel_bnb_lowE_output_tree->Branch("y_reco",&reco_var[24]);
  sel_bnb_lowE_output_tree->Branch("z_reco",&reco_var[25]);
  sel_bnb_lowE_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_bnb_lowE_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_bnb_lowE_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_bnb_lowE_output_tree->Branch("shr_charge_ratio_reco",&reco_var[29]);
  sel_bnb_lowE_output_tree->Branch("shr_consistency_reco",&reco_var[30]);
  sel_bnb_lowE_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[31]);
  sel_bnb_lowE_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[32]);  
  sel_bnb_lowE_output_tree->Branch("proton_length_reco",&reco_var[33]);
  sel_bnb_lowE_output_tree->Branch("lepton_length_reco",&reco_var[34]);
  sel_bnb_lowE_output_tree->Branch("dllee_scedr",&reco_var[38]);
  sel_bnb_lowE_output_tree->Branch("dllee_pot_weight",&reco_var[39]);
  sel_bnb_lowE_output_tree->Branch("dllee_pi0_flag_new",&reco_var[42]);
  sel_bnb_lowE_output_tree->Branch("dllee_pi0_flag_old",&reco_var[43]);
  sel_bnb_lowE_output_tree->Branch("dllee_cut_level",&reco_var[45]);
  sel_bnb_lowE_output_tree->Branch("dllee_nBDTs",&reco_var[46]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[47]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score00",&reco_var[48]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[49]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score01",&reco_var[50]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[51]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score02",&reco_var[52]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[53]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score03",&reco_var[54]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[55]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score04",&reco_var[56]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[57]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score05",&reco_var[58]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[59]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score06",&reco_var[60]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[61]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score07",&reco_var[62]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[63]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score08",&reco_var[64]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[65]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score09",&reco_var[66]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight10",&reco_var[67]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score10",&reco_var[68]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight11",&reco_var[69]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score11",&reco_var[70]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight12",&reco_var[71]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score12",&reco_var[72]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight13",&reco_var[73]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score13",&reco_var[74]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight14",&reco_var[75]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score14",&reco_var[76]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight15",&reco_var[72]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score15",&reco_var[78]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight16",&reco_var[79]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score16",&reco_var[80]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight17",&reco_var[81]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score17",&reco_var[82]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight18",&reco_var[83]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score18",&reco_var[74]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_tvweight19",&reco_var[85]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score19",&reco_var[86]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score_avg",&reco_var[87]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score_median",&reco_var[88]);
  sel_bnb_lowE_output_tree->Branch("dllee_bdt_score_max",&reco_var[89]);
  sel_bnb_lowE_output_tree->Branch("dllee_pi0_weight",&reco_var[92]);
  // ... for intrinsic nue
  TTree* sel_nue_output_tree = new TTree("sel_nue_tree", "events selected from intrinsic nue overlay sample");
  sel_nue_output_tree->Branch("run",&run);
  sel_nue_output_tree->Branch("subrun",&subrun);
  sel_nue_output_tree->Branch("event",&event);
  sel_nue_output_tree->Branch("nu_pdg",&nu_pdg);
  sel_nue_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel_nue_output_tree->Branch("nu_L_true",&nu_L_true);
  sel_nue_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel_nue_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel_nue_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel_nue_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel_nue_output_tree->Branch("tune1_weight",&tune1_weight);
  sel_nue_output_tree->Branch("spline_weight",&spline_weight);
  sel_nue_output_tree->Branch("rootino_weight",&rootino_weight);
  sel_nue_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel_nue_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel_nue_output_tree->Branch("lee_weight",&lee_weight);
  sel_nue_output_tree->Branch("sys_weights",&sys_weights);
  sel_nue_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_nue_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_nue_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_nue_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_nue_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_nue_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_nue_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_nue_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_nue_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_nue_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_nue_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_nue_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_nue_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_nue_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_nue_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_nue_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_nue_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_nue_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_nue_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_nue_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_nue_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_nue_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_nue_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_nue_output_tree->Branch("x_reco",&reco_var[23]);
  sel_nue_output_tree->Branch("y_reco",&reco_var[24]);
  sel_nue_output_tree->Branch("z_reco",&reco_var[25]);
  sel_nue_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_nue_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_nue_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_nue_output_tree->Branch("shr_charge_ratio_reco",&reco_var[29]);
  sel_nue_output_tree->Branch("shr_consistency_reco",&reco_var[30]);
  sel_nue_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[31]);
  sel_nue_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[32]);
  sel_nue_output_tree->Branch("proton_length_reco",&reco_var[33]);
  sel_nue_output_tree->Branch("lepton_length_reco",&reco_var[34]);
  sel_nue_output_tree->Branch("dllee_scedr",&reco_var[38]);
  sel_nue_output_tree->Branch("dllee_pot_weight",&reco_var[39]);
  sel_nue_output_tree->Branch("dllee_pi0_flag_new",&reco_var[42]);
  sel_nue_output_tree->Branch("dllee_pi0_flag_old",&reco_var[43]);
  sel_nue_output_tree->Branch("dllee_cut_level",&reco_var[45]);
  sel_nue_output_tree->Branch("dllee_nBDTs",&reco_var[46]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[47]);
  sel_nue_output_tree->Branch("dllee_bdt_score00",&reco_var[48]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[49]);
  sel_nue_output_tree->Branch("dllee_bdt_score01",&reco_var[50]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[51]);
  sel_nue_output_tree->Branch("dllee_bdt_score02",&reco_var[52]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[53]);
  sel_nue_output_tree->Branch("dllee_bdt_score03",&reco_var[54]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[55]);
  sel_nue_output_tree->Branch("dllee_bdt_score04",&reco_var[56]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[57]);
  sel_nue_output_tree->Branch("dllee_bdt_score05",&reco_var[58]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[59]);
  sel_nue_output_tree->Branch("dllee_bdt_score06",&reco_var[60]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[61]);
  sel_nue_output_tree->Branch("dllee_bdt_score07",&reco_var[62]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[63]);
  sel_nue_output_tree->Branch("dllee_bdt_score08",&reco_var[64]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[65]);
  sel_nue_output_tree->Branch("dllee_bdt_score09",&reco_var[66]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight10",&reco_var[67]);
  sel_nue_output_tree->Branch("dllee_bdt_score10",&reco_var[68]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight11",&reco_var[69]);
  sel_nue_output_tree->Branch("dllee_bdt_score11",&reco_var[70]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight12",&reco_var[71]);
  sel_nue_output_tree->Branch("dllee_bdt_score12",&reco_var[72]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight13",&reco_var[73]);
  sel_nue_output_tree->Branch("dllee_bdt_score13",&reco_var[74]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight14",&reco_var[75]);
  sel_nue_output_tree->Branch("dllee_bdt_score14",&reco_var[76]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight15",&reco_var[77]);
  sel_nue_output_tree->Branch("dllee_bdt_score15",&reco_var[78]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight16",&reco_var[79]);
  sel_nue_output_tree->Branch("dllee_bdt_score16",&reco_var[80]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight17",&reco_var[81]);
  sel_nue_output_tree->Branch("dllee_bdt_score17",&reco_var[82]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight18",&reco_var[83]);
  sel_nue_output_tree->Branch("dllee_bdt_score18",&reco_var[84]);
  sel_nue_output_tree->Branch("dllee_bdt_tvweight19",&reco_var[85]);
  sel_nue_output_tree->Branch("dllee_bdt_score19",&reco_var[86]);
  sel_nue_output_tree->Branch("dllee_bdt_score_avg",&reco_var[87]);
  sel_nue_output_tree->Branch("dllee_bdt_score_median",&reco_var[88]);
  sel_nue_output_tree->Branch("dllee_bdt_score_max",&reco_var[89]);
  sel_nue_output_tree->Branch("dllee_pi0_weight",&reco_var[92]);
  // ... for intrinsic nue lowE
  TTree* sel_nue_lowE_output_tree = new TTree("sel_nue_lowE_tree", "events selected from intrinsic nue lowE overlay sample");
  sel_nue_lowE_output_tree->Branch("run",&run);
  sel_nue_lowE_output_tree->Branch("subrun",&subrun);
  sel_nue_lowE_output_tree->Branch("event",&event);
  sel_nue_lowE_output_tree->Branch("nu_pdg",&nu_pdg);
  sel_nue_lowE_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel_nue_lowE_output_tree->Branch("nu_L_true",&nu_L_true);
  sel_nue_lowE_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel_nue_lowE_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel_nue_lowE_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel_nue_lowE_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel_nue_lowE_output_tree->Branch("tune1_weight",&tune1_weight);
  sel_nue_lowE_output_tree->Branch("spline_weight",&spline_weight);
  sel_nue_lowE_output_tree->Branch("rootino_weight",&rootino_weight);
  sel_nue_lowE_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel_nue_lowE_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel_nue_lowE_output_tree->Branch("lee_weight",&lee_weight);
  sel_nue_lowE_output_tree->Branch("sys_weights",&sys_weights);
  sel_nue_lowE_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_nue_lowE_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_nue_lowE_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_nue_lowE_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_nue_lowE_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_nue_lowE_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_nue_lowE_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_nue_lowE_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_nue_lowE_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_nue_lowE_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_nue_lowE_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_nue_lowE_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_nue_lowE_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_nue_lowE_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_nue_lowE_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_nue_lowE_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_nue_lowE_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_nue_lowE_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_nue_lowE_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_nue_lowE_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_nue_lowE_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_nue_lowE_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_nue_lowE_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_nue_lowE_output_tree->Branch("x_reco",&reco_var[23]);
  sel_nue_lowE_output_tree->Branch("y_reco",&reco_var[24]);
  sel_nue_lowE_output_tree->Branch("z_reco",&reco_var[25]);
  sel_nue_lowE_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_nue_lowE_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_nue_lowE_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_nue_lowE_output_tree->Branch("shr_charge_ratio_reco",&reco_var[29]);
  sel_nue_lowE_output_tree->Branch("shr_consistency_reco",&reco_var[30]);
  sel_nue_lowE_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[31]);
  sel_nue_lowE_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[32]);
  sel_nue_lowE_output_tree->Branch("proton_length_reco",&reco_var[33]);
  sel_nue_lowE_output_tree->Branch("lepton_length_reco",&reco_var[34]);
  sel_nue_lowE_output_tree->Branch("dllee_scedr",&reco_var[38]);
  sel_nue_lowE_output_tree->Branch("dllee_pot_weight",&reco_var[39]);
  sel_nue_lowE_output_tree->Branch("dllee_pi0_flag_new",&reco_var[42]);
  sel_nue_lowE_output_tree->Branch("dllee_pi0_flag_old",&reco_var[43]);
  sel_nue_lowE_output_tree->Branch("dllee_cut_level",&reco_var[45]);
  sel_nue_lowE_output_tree->Branch("dllee_nBDTs",&reco_var[46]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[47]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score00",&reco_var[48]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[49]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score01",&reco_var[50]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[51]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score02",&reco_var[52]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[53]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score03",&reco_var[54]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[55]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score04",&reco_var[56]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[57]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score05",&reco_var[58]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[59]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score06",&reco_var[60]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[61]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score07",&reco_var[62]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[63]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score08",&reco_var[64]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[65]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score09",&reco_var[66]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight10",&reco_var[67]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score10",&reco_var[68]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight11",&reco_var[69]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score11",&reco_var[70]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight12",&reco_var[71]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score12",&reco_var[72]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight13",&reco_var[73]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score13",&reco_var[74]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight14",&reco_var[75]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score14",&reco_var[76]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight15",&reco_var[77]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score15",&reco_var[78]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight16",&reco_var[79]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score16",&reco_var[80]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight17",&reco_var[81]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score17",&reco_var[82]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight18",&reco_var[83]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score18",&reco_var[84]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_tvweight19",&reco_var[85]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score19",&reco_var[86]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score_avg",&reco_var[87]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score_median",&reco_var[88]);
  sel_nue_lowE_output_tree->Branch("dllee_bdt_score_max",&reco_var[89]);
  sel_nue_lowE_output_tree->Branch("dllee_pi0_weight",&reco_var[92]);
  // ... for NC pi0
  TTree* sel_ncpi0_output_tree = new TTree("sel_ncpi0_tree", "events selected from NC pi0 overlay sample");
  sel_ncpi0_output_tree->Branch("run",&run);
  sel_ncpi0_output_tree->Branch("subrun",&subrun);
  sel_ncpi0_output_tree->Branch("event",&event);
  sel_ncpi0_output_tree->Branch("nu_pdg",&nu_pdg);
  sel_ncpi0_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel_ncpi0_output_tree->Branch("nu_L_true",&nu_L_true);
  sel_ncpi0_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel_ncpi0_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel_ncpi0_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel_ncpi0_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel_ncpi0_output_tree->Branch("tune1_weight",&tune1_weight);
  sel_ncpi0_output_tree->Branch("spline_weight",&spline_weight);
  sel_ncpi0_output_tree->Branch("rootino_weight",&rootino_weight);
  sel_ncpi0_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel_ncpi0_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel_ncpi0_output_tree->Branch("lee_weight",&lee_weight);
  sel_ncpi0_output_tree->Branch("sys_weights",&sys_weights);
  sel_ncpi0_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_ncpi0_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_ncpi0_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_ncpi0_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_ncpi0_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_ncpi0_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_ncpi0_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_ncpi0_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_ncpi0_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_ncpi0_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_ncpi0_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_ncpi0_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_ncpi0_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_ncpi0_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_ncpi0_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_ncpi0_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_ncpi0_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_ncpi0_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_ncpi0_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_ncpi0_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_ncpi0_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_ncpi0_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_ncpi0_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_ncpi0_output_tree->Branch("x_reco",&reco_var[23]);
  sel_ncpi0_output_tree->Branch("y_reco",&reco_var[24]);
  sel_ncpi0_output_tree->Branch("z_reco",&reco_var[25]);
  sel_ncpi0_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_ncpi0_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_ncpi0_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_ncpi0_output_tree->Branch("shr_charge_ratio_reco",&reco_var[29]);
  sel_ncpi0_output_tree->Branch("shr_consistency_reco",&reco_var[30]);
  sel_ncpi0_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[31]);
  sel_ncpi0_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[32]);
  sel_ncpi0_output_tree->Branch("proton_length_reco",&reco_var[33]);
  sel_ncpi0_output_tree->Branch("lepton_length_reco",&reco_var[34]);
  sel_ncpi0_output_tree->Branch("dllee_scedr",&reco_var[38]);
  sel_ncpi0_output_tree->Branch("dllee_pot_weight",&reco_var[39]);
  sel_ncpi0_output_tree->Branch("dllee_pi0_flag_new",&reco_var[42]);
  sel_ncpi0_output_tree->Branch("dllee_pi0_flag_old",&reco_var[43]);
  sel_ncpi0_output_tree->Branch("dllee_cut_level",&reco_var[45]);
  sel_ncpi0_output_tree->Branch("dllee_nBDTs",&reco_var[46]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[47]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score00",&reco_var[48]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[49]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score01",&reco_var[50]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[51]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score02",&reco_var[52]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[53]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score03",&reco_var[54]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[55]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score04",&reco_var[56]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[57]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score05",&reco_var[58]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[59]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score06",&reco_var[60]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[61]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score07",&reco_var[62]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[63]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score08",&reco_var[64]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[65]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score09",&reco_var[66]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight10",&reco_var[67]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score10",&reco_var[68]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight11",&reco_var[69]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score11",&reco_var[70]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight12",&reco_var[71]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score12",&reco_var[72]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight13",&reco_var[73]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score13",&reco_var[74]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight14",&reco_var[75]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score14",&reco_var[76]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight15",&reco_var[77]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score15",&reco_var[78]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight16",&reco_var[79]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score16",&reco_var[80]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight17",&reco_var[81]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score17",&reco_var[82]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight18",&reco_var[83]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score18",&reco_var[84]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_tvweight19",&reco_var[85]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score19",&reco_var[86]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score_avg",&reco_var[87]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score_median",&reco_var[88]);
  sel_ncpi0_output_tree->Branch("dllee_bdt_score_max",&reco_var[89]);
  sel_ncpi0_output_tree->Branch("dllee_pi0_weight",&reco_var[92]);
  // ... for CC pi0
  TTree* sel_ccpi0_output_tree = new TTree("sel_ccpi0_tree", "events selected from CC pi0 overlay sample");
  sel_ccpi0_output_tree->Branch("run",&run);
  sel_ccpi0_output_tree->Branch("subrun",&subrun);
  sel_ccpi0_output_tree->Branch("event",&event);
  sel_ccpi0_output_tree->Branch("nu_pdg",&nu_pdg);
  sel_ccpi0_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel_ccpi0_output_tree->Branch("nu_L_true",&nu_L_true);
  sel_ccpi0_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel_ccpi0_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel_ccpi0_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel_ccpi0_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel_ccpi0_output_tree->Branch("tune1_weight",&tune1_weight);
  sel_ccpi0_output_tree->Branch("spline_weight",&spline_weight);
  sel_ccpi0_output_tree->Branch("rootino_weight",&rootino_weight);
  sel_ccpi0_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel_ccpi0_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel_ccpi0_output_tree->Branch("lee_weight",&lee_weight);
  sel_ccpi0_output_tree->Branch("sys_weights",&sys_weights);
  sel_ccpi0_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_ccpi0_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_ccpi0_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_ccpi0_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_ccpi0_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_ccpi0_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_ccpi0_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_ccpi0_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_ccpi0_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_ccpi0_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_ccpi0_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_ccpi0_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_ccpi0_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_ccpi0_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_ccpi0_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_ccpi0_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_ccpi0_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_ccpi0_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_ccpi0_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_ccpi0_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_ccpi0_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_ccpi0_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_ccpi0_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_ccpi0_output_tree->Branch("x_reco",&reco_var[23]);
  sel_ccpi0_output_tree->Branch("y_reco",&reco_var[24]);
  sel_ccpi0_output_tree->Branch("z_reco",&reco_var[25]);
  sel_ccpi0_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_ccpi0_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_ccpi0_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_ccpi0_output_tree->Branch("shr_charge_ratio_reco",&reco_var[29]);
  sel_ccpi0_output_tree->Branch("shr_consistency_reco",&reco_var[30]);
  sel_ccpi0_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[31]);
  sel_ccpi0_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[32]);
  sel_ccpi0_output_tree->Branch("proton_length_reco",&reco_var[33]);
  sel_ccpi0_output_tree->Branch("lepton_length_reco",&reco_var[34]);
  sel_ccpi0_output_tree->Branch("dllee_scedr",&reco_var[38]);
  sel_ccpi0_output_tree->Branch("dllee_pot_weight",&reco_var[39]);
  sel_ccpi0_output_tree->Branch("dllee_pi0_flag_new",&reco_var[42]);
  sel_ccpi0_output_tree->Branch("dllee_pi0_flag_old",&reco_var[43]);
  sel_ccpi0_output_tree->Branch("dllee_cut_level",&reco_var[45]);
  sel_ccpi0_output_tree->Branch("dllee_nBDTs",&reco_var[46]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[47]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score00",&reco_var[48]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[49]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score01",&reco_var[50]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[51]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score02",&reco_var[52]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[53]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score03",&reco_var[54]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[55]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score04",&reco_var[56]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[57]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score05",&reco_var[58]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[59]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score06",&reco_var[60]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[61]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score07",&reco_var[62]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[63]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score08",&reco_var[54]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[65]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score09",&reco_var[66]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight10",&reco_var[67]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score10",&reco_var[68]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight11",&reco_var[69]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score11",&reco_var[70]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight12",&reco_var[71]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score12",&reco_var[72]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight13",&reco_var[73]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score13",&reco_var[74]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight14",&reco_var[75]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score14",&reco_var[76]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight15",&reco_var[77]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score15",&reco_var[78]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight16",&reco_var[79]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score16",&reco_var[80]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight17",&reco_var[81]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score17",&reco_var[82]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight18",&reco_var[83]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score18",&reco_var[84]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_tvweight19",&reco_var[85]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score19",&reco_var[86]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score_avg",&reco_var[87]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score_median",&reco_var[88]);
  sel_ccpi0_output_tree->Branch("dllee_bdt_score_max",&reco_var[89]);
  sel_ccpi0_output_tree->Branch("dllee_pi0_weight",&reco_var[92]);
  // ... for dirt nu
  TTree* sel_dirt_output_tree = new TTree("sel_dirt_tree", "events selected from dirt nu overlay sample");
  sel_dirt_output_tree->Branch("run",&run);
  sel_dirt_output_tree->Branch("subrun",&subrun);
  sel_dirt_output_tree->Branch("event",&event);
  sel_dirt_output_tree->Branch("nu_pdg",&nu_pdg);
  sel_dirt_output_tree->Branch("nu_energy_true",&nu_energy_true);
  sel_dirt_output_tree->Branch("nu_L_true",&nu_L_true);
  sel_dirt_output_tree->Branch("nu_interaction_ccnc",&nu_interaction_ccnc);
  sel_dirt_output_tree->Branch("nu_interaction_mode",&nu_interaction_mode);
  sel_dirt_output_tree->Branch("nu_interaction_type",&nu_interaction_type);
  sel_dirt_output_tree->Branch("nu_target_pdg",&nu_target_pdg);
  sel_dirt_output_tree->Branch("tune1_weight",&tune1_weight);
  sel_dirt_output_tree->Branch("spline_weight",&spline_weight);
  sel_dirt_output_tree->Branch("rootino_weight",&rootino_weight);
  sel_dirt_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel_dirt_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel_dirt_output_tree->Branch("lee_weight",&lee_weight);
  sel_dirt_output_tree->Branch("sys_weights",&sys_weights);
  sel_dirt_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_dirt_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_dirt_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_dirt_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_dirt_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_dirt_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_dirt_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_dirt_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_dirt_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_dirt_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_dirt_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_dirt_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_dirt_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_dirt_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_dirt_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_dirt_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_dirt_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_dirt_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_dirt_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_dirt_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_dirt_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_dirt_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_dirt_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_dirt_output_tree->Branch("x_reco",&reco_var[23]);
  sel_dirt_output_tree->Branch("y_reco",&reco_var[24]);
  sel_dirt_output_tree->Branch("z_reco",&reco_var[25]);
  sel_dirt_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_dirt_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_dirt_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_dirt_output_tree->Branch("shr_charge_ratio_reco",&reco_var[29]);
  sel_dirt_output_tree->Branch("shr_consistency_reco",&reco_var[30]);
  sel_dirt_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[31]);
  sel_dirt_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[32]);
  sel_dirt_output_tree->Branch("proton_length_reco",&reco_var[33]);
  sel_dirt_output_tree->Branch("lepton_length_reco",&reco_var[34]);
  sel_dirt_output_tree->Branch("dllee_scedr",&reco_var[38]);
  sel_dirt_output_tree->Branch("dllee_pot_weight",&reco_var[39]);
  sel_dirt_output_tree->Branch("dllee_pi0_flag_new",&reco_var[42]);
  sel_dirt_output_tree->Branch("dllee_pi0_flag_old",&reco_var[43]);
  sel_dirt_output_tree->Branch("dllee_cut_level",&reco_var[45]);
  sel_dirt_output_tree->Branch("dllee_nBDTs",&reco_var[46]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[47]);
  sel_dirt_output_tree->Branch("dllee_bdt_score00",&reco_var[48]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[49]);
  sel_dirt_output_tree->Branch("dllee_bdt_score01",&reco_var[50]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[51]);
  sel_dirt_output_tree->Branch("dllee_bdt_score02",&reco_var[52]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[43]);
  sel_dirt_output_tree->Branch("dllee_bdt_score03",&reco_var[54]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[55]);
  sel_dirt_output_tree->Branch("dllee_bdt_score04",&reco_var[56]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[57]);
  sel_dirt_output_tree->Branch("dllee_bdt_score05",&reco_var[58]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[59]);
  sel_dirt_output_tree->Branch("dllee_bdt_score06",&reco_var[60]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[61]);
  sel_dirt_output_tree->Branch("dllee_bdt_score07",&reco_var[62]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[63]);
  sel_dirt_output_tree->Branch("dllee_bdt_score08",&reco_var[64]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[65]);
  sel_dirt_output_tree->Branch("dllee_bdt_score09",&reco_var[66]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight10",&reco_var[67]);
  sel_dirt_output_tree->Branch("dllee_bdt_score10",&reco_var[68]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight11",&reco_var[69]);
  sel_dirt_output_tree->Branch("dllee_bdt_score11",&reco_var[70]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight12",&reco_var[71]);
  sel_dirt_output_tree->Branch("dllee_bdt_score12",&reco_var[72]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight13",&reco_var[73]);
  sel_dirt_output_tree->Branch("dllee_bdt_score13",&reco_var[74]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight14",&reco_var[75]);
  sel_dirt_output_tree->Branch("dllee_bdt_score14",&reco_var[76]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight15",&reco_var[77]);
  sel_dirt_output_tree->Branch("dllee_bdt_score15",&reco_var[78]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight16",&reco_var[79]);
  sel_dirt_output_tree->Branch("dllee_bdt_score16",&reco_var[80]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight17",&reco_var[81]);
  sel_dirt_output_tree->Branch("dllee_bdt_score17",&reco_var[82]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight18",&reco_var[83]);
  sel_dirt_output_tree->Branch("dllee_bdt_score18",&reco_var[84]);
  sel_dirt_output_tree->Branch("dllee_bdt_tvweight19",&reco_var[85]);
  sel_dirt_output_tree->Branch("dllee_bdt_score19",&reco_var[86]);
  sel_dirt_output_tree->Branch("dllee_bdt_score_avg",&reco_var[87]);
  sel_dirt_output_tree->Branch("dllee_bdt_score_median",&reco_var[88]);
  sel_dirt_output_tree->Branch("dllee_bdt_score_max",&reco_var[89]);
  sel_dirt_output_tree->Branch("dllee_pi0_weight",&reco_var[92]);
  // ... for extbnb
  TTree* sel_extbnb_output_tree = new TTree("sel_extbnb_tree", "events selected from EXTBNB data sample");
  sel_extbnb_output_tree->Branch("run",&run);
  sel_extbnb_output_tree->Branch("subrun",&subrun);
  sel_extbnb_output_tree->Branch("event",&event); 
  sel_extbnb_output_tree->Branch("tune1_weight",&tune1_weight);
  sel_extbnb_output_tree->Branch("spline_weight",&spline_weight);
  sel_extbnb_output_tree->Branch("rootino_weight",&rootino_weight);
  sel_extbnb_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel_extbnb_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel_extbnb_output_tree->Branch("lee_weight",&lee_weight);
  sel_extbnb_output_tree->Branch("sys_weights",&sys_weights);
  sel_extbnb_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_extbnb_output_tree->Branch("eta_reco",&reco_var[1]);
  sel_extbnb_output_tree->Branch("pT_reco",&reco_var[2]);
  sel_extbnb_output_tree->Branch("alphaT_reco",&reco_var[3]);
  sel_extbnb_output_tree->Branch("sphB_reco",&reco_var[4]);
  sel_extbnb_output_tree->Branch("pzEnu_reco",&reco_var[5]);
  sel_extbnb_output_tree->Branch("charge_near_trunk_reco",&reco_var[6]);
  sel_extbnb_output_tree->Branch("Q0_reco",&reco_var[7]);
  sel_extbnb_output_tree->Branch("Q3_reco",&reco_var[8]);
  sel_extbnb_output_tree->Branch("sum_thetas_reco",&reco_var[9]);
  sel_extbnb_output_tree->Branch("sum_phis_reco",&reco_var[10]);
  sel_extbnb_output_tree->Branch("pT_ratio_reco",&reco_var[11]);
  sel_extbnb_output_tree->Branch("proton_theta_reco",&reco_var[12]);
  sel_extbnb_output_tree->Branch("proton_phi_reco",&reco_var[13]);
  sel_extbnb_output_tree->Branch("min_shr_frac_reco",&reco_var[14]);
  sel_extbnb_output_tree->Branch("max_shr_frac_reco",&reco_var[15]);
  sel_extbnb_output_tree->Branch("BjxB_reco",&reco_var[16]);
  sel_extbnb_output_tree->Branch("BjyB_reco",&reco_var[17]);
  sel_extbnb_output_tree->Branch("proton_KE_reco",&reco_var[18]);
  sel_extbnb_output_tree->Branch("lepton_KE_reco",&reco_var[19]);
  sel_extbnb_output_tree->Branch("lepton_theta_reco",&reco_var[20]);
  sel_extbnb_output_tree->Branch("lepton_phi_reco",&reco_var[21]);
  sel_extbnb_output_tree->Branch("openang_reco",&reco_var[22]);
  sel_extbnb_output_tree->Branch("x_reco",&reco_var[23]);
  sel_extbnb_output_tree->Branch("y_reco",&reco_var[24]);
  sel_extbnb_output_tree->Branch("z_reco",&reco_var[25]);
  sel_extbnb_output_tree->Branch("mpid_muon_score",&reco_var[26]);
  sel_extbnb_output_tree->Branch("mpid_proton_score",&reco_var[27]);
  sel_extbnb_output_tree->Branch("mpid_electron_score",&reco_var[28]);
  sel_extbnb_output_tree->Branch("shr_charge_ratio_reco",&reco_var[29]);
  sel_extbnb_output_tree->Branch("shr_consistency_reco",&reco_var[30]);
  sel_extbnb_output_tree->Branch("nu_energy_QE_lepton_reco",&reco_var[31]);
  sel_extbnb_output_tree->Branch("nu_energy_QE_proton_reco",&reco_var[32]);
  sel_extbnb_output_tree->Branch("proton_length_reco",&reco_var[33]);
  sel_extbnb_output_tree->Branch("lepton_length_reco",&reco_var[34]);
  sel_extbnb_output_tree->Branch("dllee_scedr",&reco_var[38]);
  sel_extbnb_output_tree->Branch("dllee_pot_weight",&reco_var[39]);
  sel_extbnb_output_tree->Branch("dllee_pi0_flag_new",&reco_var[42]);
  sel_extbnb_output_tree->Branch("dllee_pi0_flag_old",&reco_var[43]);
  sel_extbnb_output_tree->Branch("dllee_cut_level",&reco_var[45]);
  sel_extbnb_output_tree->Branch("dllee_nBDTs",&reco_var[46]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight00",&reco_var[47]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score00",&reco_var[48]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight01",&reco_var[49]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score01",&reco_var[50]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight02",&reco_var[51]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score02",&reco_var[52]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight03",&reco_var[53]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score03",&reco_var[54]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight04",&reco_var[55]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score04",&reco_var[56]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight05",&reco_var[57]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score05",&reco_var[58]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight06",&reco_var[59]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score06",&reco_var[60]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight07",&reco_var[61]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score07",&reco_var[62]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight08",&reco_var[63]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score08",&reco_var[64]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight09",&reco_var[65]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score09",&reco_var[66]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight10",&reco_var[67]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score10",&reco_var[68]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight11",&reco_var[69]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score11",&reco_var[70]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight12",&reco_var[71]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score12",&reco_var[72]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight13",&reco_var[73]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score13",&reco_var[74]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight14",&reco_var[75]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score14",&reco_var[76]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight15",&reco_var[77]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score15",&reco_var[78]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight16",&reco_var[79]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score16",&reco_var[80]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight17",&reco_var[81]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score17",&reco_var[82]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight18",&reco_var[83]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score18",&reco_var[84]);
  sel_extbnb_output_tree->Branch("dllee_bdt_tvweight19",&reco_var[85]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score19",&reco_var[86]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score_avg",&reco_var[87]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score_median",&reco_var[88]);
  sel_extbnb_output_tree->Branch("dllee_bdt_score_max",&reco_var[89]);
  sel_extbnb_output_tree->Branch("dllee_pi0_weight",&reco_var[92]);
  // ... for bnb data
  TTree* sel_data_output_tree = new TTree("sel_data_tree", "events selected from BNB data sample");
  sel_data_output_tree->Branch("run",&run);
  sel_data_output_tree->Branch("subrun",&subrun);
  sel_data_output_tree->Branch("event",&event); 
  sel_data_output_tree->Branch("nu_energy_reco",&nu_energy_reco);
  sel_data_output_tree->Branch("tune1_weight",&tune1_weight);
  sel_data_output_tree->Branch("spline_weight",&spline_weight);
  sel_data_output_tree->Branch("rootino_weight",&rootino_weight);
  sel_data_output_tree->Branch("ub_tune_weight",&ub_tune_weight);
  sel_data_output_tree->Branch("xsec_corr_weight",&xsec_corr_weight);
  sel_data_output_tree->Branch("lee_weight",&lee_weight);
  sel_data_output_tree->Branch("sys_weights",&sys_weights);

  std::cout << "Initialized trees in output files" << std::endl;

  
  // * Get splines for Tune 1 reweighting * //
  
  // Open the v3.0.4 file, get the graphs
  TFile* v304_xsec_spline = new TFile(v304_xsec_spline_file, "READ");
  TGraph* xsec_mcc9_graph_numu = (TGraph*)v304_xsec_spline->Get("nu_mu_Ar40/qel_cc_n");
  TGraph* xsec_mcc9_graph_numubar = (TGraph*)v304_xsec_spline->Get("nu_mu_bar_Ar40/qel_cc_p");
  TGraph* xsec_mcc9_graph_nue = (TGraph*)v304_xsec_spline->Get("nu_e_Ar40/qel_cc_n");
  TGraph* xsec_mcc9_graph_nuebar = (TGraph*)v304_xsec_spline->Get("nu_e_bar_Ar40/qel_cc_p");
  // Open the v3.0.4a file, get the graphs
  TFile* v304a_xsec_spline = new TFile(v304a_xsec_spline_file, "READ");
  //TGraph* xsec_mcc9fixed_graph_numu = (TGraph*)v304a_xsec_spline->Get("nu_mu_Ar40/qel_cc_n");
  //TGraph* xsec_mcc9fixed_graph_numubar = (TGraph*)v304a_xsec_spline->Get("nu_mu_bar_Ar40/qel_cc_p");
  TGraph* xsec_mcc9fixed_graph_nue = (TGraph*)v304a_xsec_spline->Get("nu_e_Ar40/qel_cc_n");
  TGraph* xsec_mcc9fixed_graph_nuebar = (TGraph*)v304a_xsec_spline->Get("nu_e_bar_Ar40/qel_cc_p");
  // Open the Tune 1 file, get the graphs
  TFile* tune1_xsec_spline = new TFile(tune1_xsec_spline_file, "READ");
  TGraph* xsec_tune1_graph_numu = (TGraph*)tune1_xsec_spline->Get("nu_mu_Ar40/qel_cc_n");
  TGraph* xsec_tune1_graph_numubar = (TGraph*)tune1_xsec_spline->Get("nu_mu_bar_Ar40/qel_cc_p");
  TGraph* xsec_tune1_graph_nue = (TGraph*)tune1_xsec_spline->Get("nu_e_Ar40/qel_cc_n");
  TGraph* xsec_tune1_graph_nuebar = (TGraph*)tune1_xsec_spline->Get("nu_e_bar_Ar40/qel_cc_p");
  
  std::cout << "Got splines for Tune 1 reweighting" << std::endl;

  
  // * Fill trees in the output file * //  

  // Read in bnb nu arborist file
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
  
  // Loop over the bnb nu tree
  for ( int i=0; i < bnb_evweight_tree->GetEntries(); i++ ) {
    
    // Get the entry
    bnb_evweight_tree->GetEntry(i);
    
    // Get the rse
    rse[0] = run;
    rse[1] = subrun;
    rse[2] = event;
    
    // If there's a corresponding entry in the selected map...
    if ( sel_bnb_event_map.find(rse) != sel_bnb_event_map.end() ) {
      // Get the reconstructed variables
      reco_var = sel_bnb_event_map.find(rse)->second;
      // Check that the true energies agree... If the true energy is off by more than 0.01 MeV, this is not the right entry in the tree and we skip it
      if ( std::abs(nu_energy_true - reco_var[35]) > 0.01 ) continue;
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
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_bnb has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_bnb has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the tree
      sel_bnb_output_tree->Fill();
    }
        
  }
  
  std::cout << "Filled bnb nu tree" << std::endl;

  /*
  // Read in bnb nu lowE arborist file
  TFile* bnb_lowE_arborist = new TFile(bnb_lowE_arborist_file, "READ");
  if ( !bnb_lowE_arborist->IsOpen() ) std::cout << "Failed to open file " << bnb_lowE_arborist_file << std::endl;
  
  // Get the desired tree
  TTree* bnb_lowE_evweight_tree = (TTree*)bnb_lowE_arborist->Get("arborist/eventweight_tree");
  bnb_lowE_evweight_tree->SetBranchAddress("run",&run);
  bnb_lowE_evweight_tree->SetBranchAddress("subrun",&subrun);
  bnb_lowE_evweight_tree->SetBranchAddress("event",&event);
  bnb_lowE_evweight_tree->SetBranchAddress("nu_pdg",&nu_pdg);
  bnb_lowE_evweight_tree->SetBranchAddress("nu_energy_true",&nu_energy_true);
  bnb_lowE_evweight_tree->SetBranchAddress("nu_interaction_ccnc",&nu_interaction_ccnc);
  bnb_lowE_evweight_tree->SetBranchAddress("nu_interaction_mode",&nu_interaction_mode);
  bnb_lowE_evweight_tree->SetBranchAddress("nu_interaction_type",&nu_interaction_type);
  bnb_lowE_evweight_tree->SetBranchAddress("nu_target_pdg",&nu_target_pdg);
  bnb_lowE_evweight_tree->SetBranchAddress("nu_L_true",&nu_L_true);
  bnb_lowE_evweight_tree->SetBranchAddress("spline_weight",&spline_weight);
  bnb_lowE_evweight_tree->SetBranchAddress("rootino_weight",&rootino_weight);
  bnb_lowE_evweight_tree->SetBranchAddress("ub_tune_weight",&ub_tune_weight);
  bnb_lowE_evweight_tree->SetBranchAddress("xsec_corr_weight",&xsec_corr_weight);
  bnb_lowE_evweight_tree->SetBranchAddress("lee_weight",&lee_weight);
  bnb_lowE_evweight_tree->SetBranchAddress("sys_weights",&sys_weights_ptr);
  
  // Loop over the bnb nu lowE tree
  for ( int i=0; i < bnb_lowE_evweight_tree->GetEntries(); i++ ) {
    
    // Get the entry
    bnb_lowE_evweight_tree->GetEntry(i);
    
    // Get the rse
    rse[0] = run;
    rse[1] = subrun;
    rse[2] = event;
    
    // If there's a corresponding entry in the selected map...
    if ( sel_bnb_lowE_event_map.find(rse) != sel_bnb_lowE_event_map.end() ) {
      // Get the reconstructed energy
      reco_var = sel_bnb_lowE_event_map.find(rse)->second;
      // Get the Tune 1 weight
      if ( nu_interaction_type == 1001 ) {
	if ( nu_pdg == 14 ) tune1_weight = xsec_tune1_graph_numu->Eval(nu_energy_true/1000.) / xsec_mcc9fixed_graph_numu->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -14 ) tune1_weight = xsec_tune1_graph_numubar->Eval(nu_energy_true/1000.) / xsec_mcc9fixed_graph_numubar->Eval(nu_energy_true/1000.);
	if ( nu_pdg == 12 ) tune1_weight = xsec_tune1_graph_nue->Eval(nu_energy_true/1000.) / xsec_mcc9fixed_graph_nue->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -12 ) tune1_weight = xsec_tune1_graph_nuebar->Eval(nu_energy_true/1000.) / xsec_mcc9fixed_graph_nuebar->Eval(nu_energy_true/1000.);
      }
      else tune1_weight = 1.;
      // If the ub_tune_weight is infinite or negative, set it to 1 -- per recommendation from Steven G.
      if ( std::isinf(ub_tune_weight) || ub_tune_weight < 0. ) {
	if ( std::isinf(ub_tune_weight) )
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_bnb_lowE has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_bnb_lowE has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the tree
      sel_bnb_lowE_output_tree->Fill();
    }
    
  }
  
  std::cout << "Filled bnb nu lowE tree" << std::endl;
  */

  
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
    
    // If there's a corresponding entry in the selected map...
    if ( sel_nue_event_map.find(rse) != sel_nue_event_map.end() ) {
      // Get the reconstructed energy
      reco_var = sel_nue_event_map.find(rse)->second;
      // Get the Tune 1 weight
      if ( nu_interaction_type == 1001 ) {
	if ( nu_pdg == 12 ) tune1_weight = xsec_tune1_graph_nue->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nue->Eval(nu_energy_true/1000.);
	if ( nu_pdg == -12 ) tune1_weight = xsec_tune1_graph_nuebar->Eval(nu_energy_true/1000.) / xsec_mcc9_graph_nuebar->Eval(nu_energy_true/1000.);
      }
      else tune1_weight = 1.;
      // If the ub_tune_weight is infinite or negative, set it to 1 -- per recommendation from Steven G.
      if ( std::isinf(ub_tune_weight) || ub_tune_weight < 0. ) {
	if ( std::isinf(ub_tune_weight) )
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_nue has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_nue has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the output tree
      sel_nue_output_tree->Fill();
    }
    
  }
  
  std::cout << "Filled nue tree" << std::endl;

  
  // If we have a nue lowE event list...
  if ( sel_nue_lowE_event_map.size() > 0 ) {
    // Read in nue lowE arborist file
    TFile* nue_lowE_arborist = new TFile(nue_lowE_arborist_file, "READ");
    if ( !nue_lowE_arborist->IsOpen() ) std::cout << "Failed to open file " << nue_lowE_arborist_file << std::endl;
    
    // Get the desired tree
    TTree* nue_lowE_evweight_tree = (TTree*)nue_lowE_arborist->Get("arborist/eventweight_tree");
    // Get the desired branches
    nue_lowE_evweight_tree->SetBranchAddress("run",&run);
    nue_lowE_evweight_tree->SetBranchAddress("subrun",&subrun);
    nue_lowE_evweight_tree->SetBranchAddress("event",&event);
    nue_lowE_evweight_tree->SetBranchAddress("nu_pdg",&nu_pdg);
    nue_lowE_evweight_tree->SetBranchAddress("nu_energy_true",&nu_energy_true);
    nue_lowE_evweight_tree->SetBranchAddress("nu_interaction_ccnc",&nu_interaction_ccnc);
    nue_lowE_evweight_tree->SetBranchAddress("nu_interaction_mode",&nu_interaction_mode);
    nue_lowE_evweight_tree->SetBranchAddress("nu_interaction_type",&nu_interaction_type);
    nue_lowE_evweight_tree->SetBranchAddress("nu_target_pdg",&nu_target_pdg);
    nue_lowE_evweight_tree->SetBranchAddress("nu_L_true",&nu_L_true);
    nue_lowE_evweight_tree->SetBranchAddress("spline_weight",&spline_weight);
    nue_lowE_evweight_tree->SetBranchAddress("rootino_weight",&rootino_weight);
    nue_lowE_evweight_tree->SetBranchAddress("ub_tune_weight",&ub_tune_weight);
    nue_lowE_evweight_tree->SetBranchAddress("xsec_corr_weight",&xsec_corr_weight);
    nue_lowE_evweight_tree->SetBranchAddress("lee_weight",&lee_weight);
    nue_lowE_evweight_tree->SetBranchAddress("sys_weights",&sys_weights_ptr);
    
    // Loop over the nue lowE tree
    for ( int i=0; i < nue_lowE_evweight_tree->GetEntries(); i++ ) {
      
      // Get the entry
      nue_lowE_evweight_tree->GetEntry(i);
      
      // Get the rse
      rse[0] = run;
      rse[1] = subrun;
      rse[2] = event;
      
      // If there's a corresponding entry in the selected map...
      if ( sel_nue_lowE_event_map.find(rse) != sel_nue_lowE_event_map.end() ) {
	// Get the reconstructed energy
	reco_var = sel_nue_lowE_event_map.find(rse)->second;
	// Get the Tune 1 weight
	if ( nu_interaction_type == 1001 ) {
	  if ( nu_pdg == 12 ) tune1_weight = xsec_tune1_graph_nue->Eval(nu_energy_true/1000.) / xsec_mcc9fixed_graph_nue->Eval(nu_energy_true/1000.);
	  if ( nu_pdg == -12 ) tune1_weight = xsec_tune1_graph_nuebar->Eval(nu_energy_true/1000.) / xsec_mcc9fixed_graph_nuebar->Eval(nu_energy_true/1000.);
	}
	else tune1_weight = 1.;
	// If the ub_tune_weight is infinite or negative, set it to 1 -- per recommendation from Steven G.
	if ( std::isinf(ub_tune_weight) || ub_tune_weight < 0. ) {
	  if ( std::isinf(ub_tune_weight) )
	    std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_nue_lowE has infinite ub_tune_weight... Setting to 1" << std::endl;
	  if ( ub_tune_weight < 0. ) 
	    std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_nue_lowE has negative ub_tune_weight... Setting to 1" << std::endl;
	  ub_tune_weight = 1.;
	  xsec_corr_weight = spline_weight * rootino_weight;
	  for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
	}
	// Fill the output tree
	sel_nue_lowE_output_tree->Fill();
      }
      
    }
  
    std::cout << "Filled nue lowE tree" << std::endl;
  } // end filling nue lowE tree

  // If we have an NC pi0 event list...
  if ( sel_ncpi0_event_map.size() > 0 ) {
    // Read in NC pi0 arborist file
    TFile* ncpi0_arborist = new TFile(ncpi0_arborist_file, "READ");
    if ( !ncpi0_arborist->IsOpen() ) std::cout << "Failed to open file " << ncpi0_arborist_file << std::endl;
    
    // Get the desired tree
    TTree* ncpi0_evweight_tree = (TTree*)ncpi0_arborist->Get("arborist/eventweight_tree");
    // Get the desired branches
    ncpi0_evweight_tree->SetBranchAddress("run",&run);
    ncpi0_evweight_tree->SetBranchAddress("subrun",&subrun);
    ncpi0_evweight_tree->SetBranchAddress("event",&event);
    ncpi0_evweight_tree->SetBranchAddress("nu_pdg",&nu_pdg);
    ncpi0_evweight_tree->SetBranchAddress("nu_energy_true",&nu_energy_true);
    ncpi0_evweight_tree->SetBranchAddress("nu_interaction_ccnc",&nu_interaction_ccnc);
    ncpi0_evweight_tree->SetBranchAddress("nu_interaction_mode",&nu_interaction_mode);
    ncpi0_evweight_tree->SetBranchAddress("nu_interaction_type",&nu_interaction_type);
    ncpi0_evweight_tree->SetBranchAddress("nu_target_pdg",&nu_target_pdg);
    ncpi0_evweight_tree->SetBranchAddress("nu_L_true",&nu_L_true);
    ncpi0_evweight_tree->SetBranchAddress("spline_weight",&spline_weight);
    ncpi0_evweight_tree->SetBranchAddress("rootino_weight",&rootino_weight);
    ncpi0_evweight_tree->SetBranchAddress("ub_tune_weight",&ub_tune_weight);
    ncpi0_evweight_tree->SetBranchAddress("xsec_corr_weight",&xsec_corr_weight);
    ncpi0_evweight_tree->SetBranchAddress("lee_weight",&lee_weight);
    ncpi0_evweight_tree->SetBranchAddress("sys_weights",&sys_weights_ptr);
    
    // Loop over the NC pi0 tree
    for ( int i=0; i < ncpi0_evweight_tree->GetEntries(); i++ ) {
      
      // Get the entry
      ncpi0_evweight_tree->GetEntry(i);
      
      // Get the rse
      rse[0] = run;
      rse[1] = subrun;
      rse[2] = event;
      
      // If there's a corresponding entry in the selected map...
      if ( sel_ncpi0_event_map.find(rse) != sel_ncpi0_event_map.end() ) {
	// Get the reconstructed energy
	reco_var = sel_ncpi0_event_map.find(rse)->second;
	// Get the Tune 1 weight -- should only be NC in this sample, but just in case...
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
	    std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_ncpi0 has infinite ub_tune_weight... Setting to 1" << std::endl;
	  if ( ub_tune_weight < 0. ) 
	    std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_ncpi0 has negative ub_tune_weight... Setting to 1" << std::endl;
	  ub_tune_weight = 1.;
	  xsec_corr_weight = spline_weight * rootino_weight;
	  for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
	}
	// Fill the output tree
	sel_ncpi0_output_tree->Fill();
      }
      
    }
  
    std::cout << "Filled NC pi0 tree" << std::endl;
  } // end filling NC pi0 tree

  // If we have a CC pi0 event list...
  if ( sel_ccpi0_event_map.size() > 0 ) {
    // Read in CC pi0 arborist file
    TFile* ccpi0_arborist = new TFile(ccpi0_arborist_file, "READ");
    if ( !ccpi0_arborist->IsOpen() ) std::cout << "Failed to open file " << ccpi0_arborist_file << std::endl;
    
    // Get the desired tree
    TTree* ccpi0_evweight_tree = (TTree*)ccpi0_arborist->Get("arborist/eventweight_tree");
    // Get the desired branches
    ccpi0_evweight_tree->SetBranchAddress("run",&run);
    ccpi0_evweight_tree->SetBranchAddress("subrun",&subrun);
    ccpi0_evweight_tree->SetBranchAddress("event",&event);
    ccpi0_evweight_tree->SetBranchAddress("nu_pdg",&nu_pdg);
    ccpi0_evweight_tree->SetBranchAddress("nu_energy_true",&nu_energy_true);
    ccpi0_evweight_tree->SetBranchAddress("nu_interaction_ccnc",&nu_interaction_ccnc);
    ccpi0_evweight_tree->SetBranchAddress("nu_interaction_mode",&nu_interaction_mode);
    ccpi0_evweight_tree->SetBranchAddress("nu_interaction_type",&nu_interaction_type);
    ccpi0_evweight_tree->SetBranchAddress("nu_target_pdg",&nu_target_pdg);
    ccpi0_evweight_tree->SetBranchAddress("nu_L_true",&nu_L_true);
    ccpi0_evweight_tree->SetBranchAddress("spline_weight",&spline_weight);
    ccpi0_evweight_tree->SetBranchAddress("rootino_weight",&rootino_weight);
    ccpi0_evweight_tree->SetBranchAddress("ub_tune_weight",&ub_tune_weight);
    ccpi0_evweight_tree->SetBranchAddress("xsec_corr_weight",&xsec_corr_weight);
    ccpi0_evweight_tree->SetBranchAddress("lee_weight",&lee_weight);
    ccpi0_evweight_tree->SetBranchAddress("sys_weights",&sys_weights_ptr);
    
    // Loop over the CC pi0 tree
    for ( int i=0; i < ccpi0_evweight_tree->GetEntries(); i++ ) {
      
      // Get the entry
      ccpi0_evweight_tree->GetEntry(i);
      
      // Get the rse
      rse[0] = run;
      rse[1] = subrun;
      rse[2] = event;
      
      // If there's a corresponding entry in the selected map...
      if ( sel_ccpi0_event_map.find(rse) != sel_ccpi0_event_map.end() ) {
	// Get the reconstructed energy
	reco_var = sel_ccpi0_event_map.find(rse)->second;
	// Get the Tune 1 weight -- should only be nu_pdg==14 in this sample, but just do everything...
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
	    std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_ccpi0 has infinite ub_tune_weight... Setting to 1" << std::endl;
	  if ( ub_tune_weight < 0. ) 
	    std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_ccpi0 has negative ub_tune_weight... Setting to 1" << std::endl;
	  ub_tune_weight = 1.;
	  xsec_corr_weight = spline_weight * rootino_weight;
	  for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
	}
	// Fill the output tree
	sel_ccpi0_output_tree->Fill();
      }
      
    }
  
    std::cout << "Filled CC pi0 tree" << std::endl;
  } // end filling CC pi0 tree

  /*
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
    
    // If there's a corresponding entry in the selected map...
    if ( sel_dirt_event_map.find(rse) != sel_dirt_event_map.end() ) {
      // Get the reconstructed energy
      reco_var = sel_dirt_event_map.find(rse)->second;
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
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_dirt has infinite ub_tune_weight... Setting to 1" << std::endl;
	if ( ub_tune_weight < 0. ) 
	  std::cout << "Event (" << run << ", " << subrun << ", " << event << ") in sel_dirt has negative ub_tune_weight... Setting to 1" << std::endl;
	ub_tune_weight = 1.;
	xsec_corr_weight = spline_weight * rootino_weight;
	for (const auto& sys_weight_pair : sys_weights) { for (auto weight : sys_weight_pair.second) { weight = 1.; } }
      }
      // Fill the tree
      sel_dirt_output_tree->Fill();
    }    
    
  }
  
  std::cout << "Filled dirt tree" << std::endl;
  */

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
  
  // Fill extbnb tree
  for ( auto event_it = sel_extbnb_event_map.begin(); event_it != sel_extbnb_event_map.end(); event_it++ ) {
    // Get the variables to be put into the tree                                                                                                                                 
    run = event_it->first[0];
    subrun = event_it->first[1];
    event = event_it->first[2];
    reco_var = event_it->second;
    // Fill the tree                                                                                                                                  
    sel_extbnb_output_tree->Fill();
  }
  
  // Fill bnb data tree
  for ( auto event_it = sel_data_event_map.begin(); event_it != sel_data_event_map.end(); event_it++ ) {
    // Get the variables to be put into the tree                                                                                                                                 
    run = event_it->first[0];
    subrun = event_it->first[1];
    event = event_it->first[2];
    nu_energy_reco = event_it->second;
    // Fill the tree                                                                                                                                                             
    sel_data_output_tree->Fill();
  }

  std::cout << "Filled data trees" << std::endl;


  // * Conclusion * //
  
  // Count events in all output trees and make sure we did this right
  int n_sel_bnb_in       = sel_bnb_event_map.size();
  int n_sel_bnb_out      = sel_bnb_output_tree->GetEntries();
  int n_sel_bnb_lowE_in  = sel_bnb_lowE_event_map.size();
  int n_sel_bnb_lowE_out = sel_bnb_lowE_output_tree->GetEntries();
  int n_sel_nue_in       = sel_nue_event_map.size();
  int n_sel_nue_out      = sel_nue_output_tree->GetEntries();
  int n_sel_nue_lowE_in  = sel_nue_lowE_event_map.size();
  int n_sel_nue_lowE_out = sel_nue_lowE_output_tree->GetEntries();
  int n_sel_ncpi0_in     = sel_ncpi0_event_map.size();
  int n_sel_ncpi0_out    = sel_ncpi0_output_tree->GetEntries();
  int n_sel_ccpi0_in     = sel_ccpi0_event_map.size();
  int n_sel_ccpi0_out    = sel_ccpi0_output_tree->GetEntries();
  int n_sel_dirt_in      = sel_dirt_event_map.size();
  int n_sel_dirt_out     = sel_dirt_output_tree->GetEntries();
  int n_sel_extbnb_in    = sel_extbnb_event_map.size();
  int n_sel_extbnb_out   = sel_extbnb_output_tree->GetEntries();
  int n_sel_data_in      = sel_data_event_map.size();
  int n_sel_data_out     = sel_data_output_tree->GetEntries();
  
  // Print information
  // ... for MC
  std::cout << "Selected event count information" << std::endl;
  std::cout << "  BNB nu overlay events: " << n_sel_bnb_in << " events in, " << n_sel_bnb_out << " events out" << std::endl;
  std::cout << "  BNB nu lowE overlay events: " << n_sel_bnb_lowE_in << " events in, " << n_sel_bnb_lowE_out << " events out" << std::endl;
  std::cout << "  intrinsic nue overlay events: " << n_sel_nue_in << " events in, " << n_sel_nue_out << " events out" << std::endl;
  std::cout << "  intrinsic nue lowE overlay events: " << n_sel_nue_lowE_in << " events in, " << n_sel_nue_lowE_out << " events out" << std::endl;
  std::cout << "  NC pi0 overlay events: " << n_sel_ncpi0_in << " events in, " << n_sel_ncpi0_out << " events out" << std::endl;
  std::cout << "  CC pi0 overlay events: " << n_sel_ccpi0_in << " events in, " << n_sel_ccpi0_out << " events out" << std::endl;
  std::cout << "  dirt overlay events: " << n_sel_dirt_in << " events in, " << n_sel_dirt_out << " events out" << std::endl;
  std::cout << "  EXTBNB data events: " << n_sel_extbnb_in << " events in, " << n_sel_extbnb_out << " events out" << std::endl;
  if ( (n_sel_bnb_in == n_sel_bnb_out) and (n_sel_bnb_lowE_in == n_sel_bnb_lowE_out)
       and (n_sel_nue_in == n_sel_nue_out) and (n_sel_nue_lowE_in == n_sel_nue_lowE_out)
       and (n_sel_ncpi0_in == n_sel_ncpi0_out) and (n_sel_ccpi0_in == n_sel_ccpi0_out)
       and (n_sel_dirt_in == n_sel_dirt_out) and (n_sel_extbnb_in == n_sel_extbnb_out) ) {
    std::cout << "These values are consistent!" << std::endl;
  }
  else { std::cout << "Warning: Number of events in output trees is inconsistent with input lists" << std::endl; }
  // ... for data
  std::cout << "Selected in BNB data event count information" << std::endl;
  std::cout << "  BNB data events: " << n_sel_data_in << " events in, " << n_sel_data_out << " events out" << std::endl;
  if ( n_sel_data_in == n_sel_data_out ) { std::cout << "These values are consistent!" << std::endl; }
  else { std::cout << "Warning: Number of events in BNB data output tree is inconsistent with input list" << std::endl; }
  
  // Write and close the output file
  out->Write();
  out->Close();
  
  // El Fin
  std::cout << "Wrote output file " << output_file << std::endl;
  std::cout << "Done!" << std::endl;
  return 0;
  
}
