#ifndef SBNSINGLEPHOTON_H_
#define SBNSINGLEPHOTON_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNchi.h"
#include "SBNconfig.h"
#include "SBNgenerate.h"

#include "TH1.h"
#include "TH2.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TList.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMarker.h"
#include "TObjArray.h"

#include "TMath.h"
#include <ctime>
#include "params.h"

#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"

#include "Math/ProbFunc.h"
#include "Math/DistFunc.h"

#include "prob.h"
#include "ngrid.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>


namespace sbn{
 /***********
 *  A singlephoton class built purely for perform NCpi0 normalization fit, NC delta fit for gLEE analysis
 *  	Guanqun Ge, July 19 2020
 */

    class SBNsinglephoton: public SBNconfig{


	//regular grid, for a flat normalization
        NGrid m_grid;
        int m_num_total_gridpoints;
        std::vector<std::vector<double>> m_vec_grid;

        // extra grid for polynomial scaling purpose
	NGrid m_poly_grid; 
	int m_poly_total_gridpoints;
	bool m_bool_poly_grid;  //use polynomial grid or not
	std::vector<std::vector<double>> m_vec_poly_grid;


	int m_total_gridpoints;  //total number of grid points 

        std::vector<SBNspec> m_scaled_spec_grid;  

        TMatrixT<double> * m_full_fractional_covariance_matrix; //full systematic( include detector) covariance matrix
        TMatrixT<double> * m_full_but_genie_fractional_covariance_matrix;   //full systematic covariance matrix other than genie
        TMatrixT<double> * m_genie_fractional_covariance_matrix;   //genie only covariance matrix

	SBNspec* m_cv_spectrum;   //genie CV spectra
	SBNspec* m_data_spectrum; // data spectra
        SBNchi *m_chi;

	double m_cv_delta_scaling;   //NCdelta scaling,should only be applied to FAKE DATA for sensitivity study!
	bool m_bool_modify_cv;  //if CV is modified for systematic covar matrix calculation
	bool m_bool_cv_spectrum_generated;   //if CV spec is generated
	bool m_bool_cv_spectrum_loaded;   //if CV spec is loaded, for fit purposes
	bool m_bool_data_spectrum_loaded;  //if data spec is loaded

        bool m_use_CNP;

        int m_max_number_iterations;
        double m_chi_min_convergance_tolerance;

        double m_random_seed;
        std::string tag;   //tag for the analysis: NCpi0, NCDelta etc

	std::map<int, std::vector<double>> m_map;   // map of best-fit point index and chi2 vector

	//things used only when generating pre-scaled files
	int num_files;
        std::vector<int> nentries;
        std::vector<TFile *> files;
        std::vector<TTree *> trees;

        public:

	SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid);
	SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid, NGrid in_polygrid, bool has_polygrid);

	// MEMBER FUNCTION//
        int PreScaleSpectrum(std::string xmlname, double, std::vector<double>& param);	
        int PreScaleSpectrum(std::string xmlname, std::vector<double>& param);  //generate scaled spectrum per set of parameter
	int GeneratePreScaledSpectra();    //generate pre-scaled spectra for full polynomial grid
	int LoadSpectraApplyFullScaling();  //load spectra and apply full scaling
	//calc scale factors based on event energy and polynomial parameters
        double ScaleFactor(double E, double, std::vector<double>& param);
	//int WriteOutCV(std::string tag); 


	int CalcChiGridScanShapeOnlyFit();  //calculate chi2 surface for NCpi0 fit.
	int CalcChiGridScan();  //simply grid scan with one systematic covariance matrix

	int SetStatOnly();
	int SetFlatFullFracCovarianceMatrix(double );
	int SetFullFractionalCovarianceMatrix(std::string filename, std::string matrix_name);
	int SetGenieFractionalCovarianceMatrix(std::string filename);
	int CalcFullButGenieFractionalCovarMatrix();
	int PrintOutFitInfo(std::map<int, std::vector<double>>& , std::string tag, bool);

	int LoadCV();
	int LoadData(std::string filename);
	int SetPolyGrid(NGrid ingrid);
	
	//not finished yet, what information do we wnat to print out
	int SaveHistogram();
	int SaveHistogram(std::map<int, std::vector<double>>& );
	TH2D* Do2DInterpolation(int, std::vector<double>& x, std::vector<double>& y, std::vector<double>& value);
	std::vector<TGraph*> FindContour(TH2D*, int n);
	int RemoveNan(TMatrixT<double>*); //remove the nan's from matrix
	
	int ModifyCV(double factor);
	int ModifyCV(double, std::vector<double> param);

	SBNspec GeneratePointSpectra(int );

	protected:

	int OpenFiles();
	int CloseFiles();	
	int PoissonFluctuation(SBNspec *);
        /*//Member Functions
        
         int UpdateInverseCovarianceMatrixCNP(size_t best_grid_point, const std::vector<float> &datavec, TMatrixT<double>& inverse_collapsed, SBNchi * helper);
         int UpdateInverseCovarianceMatrix(size_t best_grid_point, TMatrixT<double>& inverse_collapsed, SBNchi * helper);
         std::vector<double> PerformIterativeGridFit(const std::vector<float> &datavec, const size_t grid_pt, const TMatrixT<double>& inverse_background_collapsed_covariance_matrix);


        int UseCNP(){m_use_CNP = true;};
        int FullFeldmanCousins();
        int PointFeldmanCousins(size_t);
        std::vector<double> GlobalScan();
        std::vector<double> GlobalScan(int);
        std::vector<double> GlobalScan2(int);
        std::vector<double> GlobalScan(SBNspec *obs);
        int RasterScan(); 
       */ 
/*
        int LoadPreOscillatedSpectrum(int);
        int LoadPreOscillatedSpectra();

        int SetRandomSeed(double);

        int GenerateBackgroundSpectrum(); 
        int SetBackgroundSpectrum(std::string filein, std::string scale_nam, double val);
        int LoadBackgroundSpectrum();
        int LoadBackgroundSpectrum(std::string);

        int CalcSBNchis();

        int SetCoreSpectrum(std::string);
        int SetFractionalCovarianceMatrix(TMatrixT<double> *);
        int SetFractionalCovarianceMatrix(std::string, std::string);
        int SetEmptyFractionalCovarianceMatrix();
        int SetNumUniverses(int);
        int SetStatOnly();

        NeutrinoModel convert3p1(std::vector<double> ingrid);
       

        std::vector<double> getConfidenceRegion(TGraph *gmin, TGraph *gmax,double val);
        
        int AddFlatDetSystematic(double percent);
        std::vector<TGraph*> LoadFCMaps(std::string filein);
        std::vector<TGraph*> MakeFCMaps(std::string filein);
        std::vector<TGraph*> MakeFCMaps(std::string filein,size_t);
        

        std::vector<double> MakeMedianMaps(std::string filein,size_t i);



        int GenerateScaledSpectra();
        std::string m_subchannel_to_scale;

        //This is a stopgap for better SBNchi integration.Hrump, need to fix that wierd float oddity. 
        float CalcChi(const std::vector<float>& data, const std::vector<double>& prediction, const TMatrixT<double> & inverse_covariance_matrix );

*/




    };



}
#endif
