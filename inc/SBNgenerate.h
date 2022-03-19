#ifndef SBNGENERATE_H_
#define SBNGENERATE_H_

#include <cmath>
#include <vector>
#include <iostream>

#include "SBNspec.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "prob.h"

#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TLine.h"

#include "TROOT.h"
#include "TRint.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THnSparse.h"
#include "TTreeFormula.h"

#include <map>
#include <ctime>
#include "params.h"

namespace sbn{

class SBNgenerate : public SBNconfig{
	
	public:
		
	SBNspec spec_central_value;	
	SBNspec* spec_cv;

	NeutrinoModel nu_model;

	SBNspec spec_osc_sin;
	SBNspec spec_osc_sinsq;
	
  SBNgenerate(std::string xmlname, NeutrinoModel inModel, bool cache_the_data=false );
  SBNgenerate(std::string xmlname, bool cache_the_data=false );
    ~SBNgenerate();

	virtual bool EventSelection(int file);
	virtual int FillHistograms(int file, int uni, double wei);
	
	int WritePrecomputedOscSpecs(std::string tag);
	int WriteCVSpec(std::string tag);
  int regenerate_osc( const NeutrinoModel& model );

	//Some checks on multisims

	//Multisim input variables
	std::vector<std::vector<double> > vars;

	int num_files;
	std::vector<TFile *> files;
	std::vector<TTree *> trees;

	std::vector<int> nentries;
	std::vector< TBranch *> * branch_weight;
        std::vector<std::map<std::string, std::vector<eweight_type> >* > f_weights;

	std::vector<std::vector<int> > vars_i;
	std::vector<std::vector<double> > vars_d;

  // CACHE VARIABLES
  // we build an event cache for each sample, associated to each entry in the files vector
  typedef struct EventCache_t {
    int num_events; ///< number of events in the sample
    int files_v_index; ///< entry in files that this cache is derived from
    int num_var_types; ///< variables to store in the cache
    std::vector< float > recovar;  ///< reco var
    std::vector< float > weight;   ///< global weight per event
    std::vector< std::vector<float> > data; ///< outer vector is variable type, inner vector is the cached value for each entry in the file.
    
  };
  bool _cache_event_data;
  std::vector< EventCache_t > event_cache_v; ///< event info cached for each sample

};


};
#endif
