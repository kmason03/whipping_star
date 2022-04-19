#include "SBNllminimizer.h"

#include <string>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

namespace sbn {

  SBNllminimizer* SBNllminimizer::_active_copy = NULL;
  
  SBNllminimizer::SBNllminimizer( std::string xml_config )
    : _osc_model(0.01,0.707,0.707), // these values are arbitrary
      _gen( xml_config, _osc_model, true ),
      nBins_e(12), // this needs to be not hard-coded
      nBins_mu(19),
      nBins(31),
      covFracSys(NULL)
  {
    
    TFile _fsys("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");
    covFracSys = (TMatrixD*)_fsys.Get("frac_covariance");
    _fsys.Close();
    
  }

  SBNllminimizer::~SBNllminimizer()
  {
    if ( covFracSys )
      delete covFracSys;
  }
  
  double SBNllminimizer::negative_likelihood_ratio( const double* par )
  {
    std::cout << _active_copy << std::endl;

    // function to calculate the chi2 between fake universe and the mc events with osc weights
    float dm4x  = par[0];
    float Ue4   = par[1];
    float Um4   = par[2];
    float e_app = 4*pow(Ue4,2)*pow(Um4,2);     // sin^2(2theta_mue)
    float e_dis = 4*pow(Ue4,2)*(1-pow(Ue4,2)); // sin^2(2theta_ee)
    float m_dis = 4*pow(Um4,2)*(1-pow(Um4,2)); // sin^2(2theta_mumu)

    // find the closest mass spectra
    // Note: this is from Katie's LL calc. We do not need to do this with cached events.
    // int lowidx, highidx;
    // double prevval = 0;
    // for(int i = 1;i<a_sinsqSpec.size();i++){
    //   // std::cout<<std::get<1>(a_sinsqSpec.at(i))<<" "<<par[0]<<" "<<prevval<<std::endl;
    //   if (par[0]==std::get<1>(a_sinsqSpec.at(i))){
    // 	lowidx = i;
    // 	highidx =i;
    // 	break;
    //   }
    //   else if (par[0]<std::get<1>(a_sinsqSpec.at(i)) && par[0] >prevval ){
    // 	lowidx = i-1;
    // 	highidx =i;
    // 	break;
    //   }
    //   else if( i == a_sinsqSpec.size()-1){
    // 	lowidx = i;
    // 	highidx =i;
    //   }
    //   else prevval = std::get<1>(a_sinsqSpec.at(i));
    // }

    // int closeidx = lowidx;
    // if (lowidx <0) lowidx =0;
    // if (lowidx<a_sinsqSpec.size()-2){
    //   double diff1 = par[0] - std::get<1>(a_sinsqSpec.at(lowidx));
    //   double diff2 = std::get<1>(a_sinsqSpec.at(highidx)) - par[0];
    //   if (diff2 <diff1) closeidx =highidx;
    // }
    // // get spectra and covar
    // SBNspec inspec = std::get<0>(a_sinsqSpec.at(closeidx));
    // SBNspec newSpec = GetOscillatedSpectra(cvSpec, inspec,e_app,e_dis, m_dis);

    // update the osc model
    _active_copy->_osc_model = NeutrinoModel( dm4x, Ue4, Um4 );
    _active_copy->_gen.regenerate_osc( _active_copy->_osc_model ); ///< regenerates spectrum with new dm4x

    // now we scale the different components
    _active_copy->_gen.spec_osc_sinsq.Scale("fullosc",0.0);
    _active_copy->_gen.spec_osc_sinsq.Scale("bnb",-1*m_dis);    
    _active_copy->_gen.spec_osc_sinsq.Scale("nue",-1*e_dis);
    _active_copy->_gen.spec_osc_sinsq.Scale("ext",0.0);
    _active_copy->_gen.spec_central_value.Scale("fullosc",0.0);    
    _active_copy->_gen.spec_central_value.Add(&(_active_copy->_gen.spec_osc_sinsq));

    TMatrixD tmp1_covFracSys( *(_active_copy->covFracSys) );
    TMatrixD tmp2_covFracSys( *(_active_copy->covFracSys) );    
    SBNchi tmpChi(_active_copy->_gen.spec_central_value, &tmp1_covFracSys );
    int b = tmpChi.FillCollapsedFractionalMatrix( &tmp2_covFracSys );
    TMatrixD tmpcov = SBNllminimizer::GetTotalCov(_active_copy->_observed_bins,
						  _active_copy->_gen.spec_central_value, tmp2_covFracSys );
    // calculate -2LLH
    double f = SBNllminimizer::GetLLHFromVector(_active_copy->_observed_bins,
						_active_copy->_gen.spec_central_value,
						tmpcov, false);
    
    return f;
  }
  
  int SBNllminimizer::doFit( std::vector<float>& obs_bins )
  {
    std::string minName =  "Minuit2";
    std::string algoName = "Scan";
    
    _active_copy = this;
    _active_copy->setObservedBinValues( obs_bins ); // we need a number of bins check
    
    ROOT::Math::Minimizer* min = 
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    ROOT::Math::Functor f(&SBNllminimizer::negative_likelihood_ratio,3); 
    min->SetFunction(f);

    delete min;
  }

  TMatrixD SBNllminimizer::GetTotalCov(const std::vector<float>& obsSpec, const SBNspec& expSpec, const TMatrixD& Mfracsys) {
    // function to take the fractional Msys and return totl Msys+Mstat
    // inputs:
    // obsSpec: "data" spectra
    // predSpec: "MC" spectra
    // Mfracsys: fractional (flux+xsec+detvar) covariance matrix
    int nBins = _active_copy->nBins;    
    TMatrixD fullcov(nBins,nBins);
    for(int i = 0; i < nBins; i++){
      for(int j = 0; j < nBins; j++){
	// first set to zero
	fullcov[i][j] = 0.0;
	// scale to the prediction
	fullcov[i][j] = (Mfracsys)[i][j]*expSpec.collapsed_vector[i]*expSpec.collapsed_vector[j];
	// add in stat errors start with CNP for "data" errors
	if(i==j){
	  if (expSpec.collapsed_vector[i] >0 ){
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

  
  float SBNllminimizer::GetLLHFromVector(const std::vector<float>& obsSpec,
					 const SBNspec& expSpec,
					 const TMatrixD& Msys,
					 bool prints)
  {
    // // function to calculate a chi2 (shape + rate)
    // // inputs:
    // // obsSpec: "data" vector
    // // predSpec: "MC" spectra
    // // Mfracsys: total (flux+xsec+detvar) covariance matrix
    float chisqTest = 0;
    int nBins = _active_copy->nBins;
    
    // inv cov for chi2calc
    TMatrixD invcov = Msys;
    invcov.Invert();
    
    // add the chi2-like part
    chisqTest = 0;
    for(int i = 0; i < nBins; i++){
      for(int j = 0; j < nBins; j++){
	// (obsi-predi)*(invcov)*(obsj-predj)
	if(i==j && prints) std::cout<<i<<" "
				    <<obsSpec[i]<<" "
				    <<expSpec.collapsed_vector[i]<<" "
				    <<Msys[i][j]<<" "
				    <<((obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]))
				    <<std::endl;
	chisqTest += (obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]);
      }
    }
    // now need ln(det(2Pi*M))
    // TMatrixD tempcov = 2*3.14159265358979323846*Msys;
    // std::cout<<"chi2: "<<chisqTest<<" det: "<<log(tempcov.Determinant())<<std::endl;
    chisqTest += log(Msys.Determinant());

    return chisqTest;
  }//end of GetLLHFromSpectra
  
}
