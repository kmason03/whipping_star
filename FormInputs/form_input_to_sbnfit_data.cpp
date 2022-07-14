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
  const char* sel_run1_eventlist_file = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_FakeData_set1_run1_eventlist.txt";
  const char* sel_run2_eventlist_file = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/empty.txt";
  const char* sel_run3_eventlist_file = "/uboone/app/users/yatesla/othersys_mcc9/input_to_sbnfit/sel1mu1p_FakeData_set1_run3_eventlist.txt";
  // ... for output
  const char* output_file = "/uboone/data/users/yatesla/othersys_mcc9/input_to_sbnfit/input_to_sbnfit_fakedata_set1_1mu1p_Apr20.root";
  // ... for flagging whether there is a vtxid between the RSE and reco energy
  const bool has_vtxid = false;

  // Initalize the number of reco variables
  const int N_var = 1;
  
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
	if ( has_vtxid ) reco_var[i] = std::stof(split_line[i+4]);  // skip vtxid after RSE
	else reco_var[i] = std::stof(split_line[i+3]);
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
  sel_run1_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_run1_output_tree->Branch("weights",&sys_weights);
  // ... for Run 2
  TTree* sel_run2_output_tree = new TTree("sel_run2_tree", "events selected from Run 2 sample");
  sel_run2_output_tree->Branch("run",&run);
  sel_run2_output_tree->Branch("subrun",&subrun);
  sel_run2_output_tree->Branch("event",&event);
  sel_run2_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_run2_output_tree->Branch("weights",&sys_weights);
  // ... for Run 3
  TTree* sel_run3_output_tree = new TTree("sel_run3_tree", "events selected from Run 3 sample");
  sel_run3_output_tree->Branch("run",&run);
  sel_run3_output_tree->Branch("subrun",&subrun);
  sel_run3_output_tree->Branch("event",&event);
  sel_run3_output_tree->Branch("nu_energy_reco",&reco_var[0]);
  sel_run3_output_tree->Branch("weights",&sys_weights);
  
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
