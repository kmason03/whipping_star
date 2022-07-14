#include "MillsFunctions.h"
// using namespace sbn;


#include <sstream>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

void testFunction(){
	std::cout << " Test Function Working Appropriately\n";
	assert (1==2);
}

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


void generate_spectra(int massid, std::string xml, std::string tag, double dm2_lowbound, double dm2_hibound, int dm2_grdpts){
	 //function to generate the different mass spectra we need
	 // prerunning this speeds up the initial grid search part and helps with plotting
	 // note: currently run twice to prevent the job being killed
	 // inputs:
	 // int massid: mass index to run
	 // outputs:
	 // none, but writes spectra in root files to current directory

	// start with a null model if mass id ==0
  if(massid==0){
    sbn::NeutrinoModel nullModel(0, 0, 0);
    sbn::SBNgenerate * bkgo = new sbn::SBNgenerate(xml,nullModel);
    sbn::SBNspec bkg = bkgo->spec_central_value;
    bkg.WriteOut(tag+"_Bkg");
  }
	double mnu;
	mnu = pow(10.,((massid+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
	//Model: mnu, ue4, um4
	//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
	//see davio's thesis eqs 4.14-16
	//we're precomputing, so we don't really care about the u's set to sin2=1
	sbn::NeutrinoModel testModel(mnu, sqrt(.5), sqrt(.5));

	// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
	sbn::SBNgenerate * gen = new sbn::SBNgenerate(xml,testModel);
	std::cout<<massid<<std::endl;
	// Write them to files
	gen->WritePrecomputedOscSpecs(tag);

	return;
}

void generate_spectra(int massid, int Um4sqid, std::string xml, std::string tag, double dm2_lowbound, double dm2_hibound, int dm2_grdpts, double Um4sq_lowbound, double Um4sq_hibound, int Um4sq_grdpts){
  // Generate SM as background
  if((massid==0) && (Um4sqid == 0)){
    sbn::NeutrinoModel nullModel(0, 0, 0);
    sbn::SBNgenerate * bkgo = new sbn::SBNgenerate(xml,nullModel);
    sbn::SBNspec bkg = bkgo->spec_central_value;
    bkg.WriteOut(tag+"_Bkg");
  }
  double mnu;
  mnu = pow(10.,((massid+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
  double Um4sq;
  Um4sq = pow(10.,((Um4sqid+.5)/float(Um4sq_grdpts)*TMath::Log10(sqrt(Um4sq_hibound)/sqrt(Um4sq_lowbound)) + TMath::Log10(sqrt(Um4sq_lowbound))));
  std::cout << "Mnu and Um4sq: "<< mnu << " " << Um4sq << "\n";
  sbn::NeutrinoModel testModel(mnu, sqrt(.5), sqrt(Um4sq));

  // on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
  bool isDecay = true;
  sbn::SBNgenerate * gen = new sbn::SBNgenerate(xml,testModel,isDecay);
  std::cout<<massid<<std::endl;
  // Write them to files
  gen->WritePrecomputedOscSpecs(tag);
}


void calcRTerm(double& detTerm, double& chiSQ, const std::vector<double>& obs, const std::vector<double>& exp, const TMatrixD& invCov, const TMatrixD& cov, int nBins){
	chiSQ = 0;
	detTerm = 0;
	double pi = 3.141592653589793238;
	for(int i = 0; i < nBins; i++){
		for(int j = 0; j < nBins; j++){
			// (obsi-predi)*(invcov)*(obsj-predj)
			chiSQ += (obs[i] - exp[i])*invCov[i][j]*(obs[j] - exp[j]);
			// std::cout << obs.collapsed_vector[i] << " " << exp.collapsed_vector[i] << " " << invCov[i][j] << " " << obs.collapsed_vector[j] << " " << exp.collapsed_vector[j] << "\n";
		}
	}
	TMatrixD covMult = 2*pi*cov;
	detTerm = log(covMult.Determinant());
	// std::cout << "detTerm, " << detTerm << ", chiSQ, " << chiSQ << "\n";
}

void scaleCovShapeRate(TMatrixD& cov, const TMatrixD& covFracCollapsed, const std::vector<double>& obs, const std::vector<double>& exp, int nBins){
	for(int i = 0; i < nBins; i++){
		for(int j = 0; j < nBins; j++){
			// remove some nans
			if(exp[i] > 0 && exp[j] >0 ){
				cov[i][j] += (covFracCollapsed)[i][j]*exp[i]*exp[j];
			}
			if(i==j){
				// cov[i][j] += thisGridCVSpec.collapsed_vector[i]; // Davio includes this
				// add in stat errors start with CNP for "data" errors
				if (exp[i] >0 && obs[i] >0){
					cov[i][j] += 3.0 / (1.0/obs[i] + 2.0/exp[i]);
					// poisson version
					// cov = (M-mu)**2 / (2*(mu - M + M*np.log(M/mu)))
					// cov[i][j] += pow(truespec_first.collapsed_vector[i]-thisGridCVSpec.collapsed_vector[i],2) / (2.0*(thisGridCVSpec.collapsed_vector[i]- truespec_first.collapsed_vector[i]+ truespec_first.collapsed_vector[i]*log(truespec_first.collapsed_vector[i]/cvSpec.collapsed_vector[i])));
				}
				else {
					cov[i][j] += exp[i]/2.0;
				}
				// mc stat error
				// cov[i][j] +=stkerr[i]*stkerr[i];
			}
		}
	}
}

void scaleCovShapeOnly(TMatrixD& cov, const TMatrixD& covFracCollapsed, const std::vector<double>& obs, std::vector<double>& exp, int nBins){
	double covIntegral(0.0), predIntegral(0.0), obsIntegral(0.0),fnorm;
	for(int i = 0; i < nBins; i++){
		predIntegral+=exp[i];
		obsIntegral +=obs[i];
		for(int j = 0; j < nBins; j++){
			covIntegral += (covFracCollapsed)[i][j]*exp[i]*exp[j];
		}
	}
  fnorm = covIntegral/pow(predIntegral,2);
	for(int i = 0; i < nBins; i++){
		exp[i] *= (obsIntegral/predIntegral); // normalize prediction
	}
	for(int i = 0; i < nBins; i++){
		for(int j = 0; j < nBins; j++){
			cov[i][j] = ((covFracCollapsed)[i][j]-fnorm)*exp[i]*exp[j];
			// std::cout << i << " " << j << " " << ((*covFracCollapsed)[i][j]-fnorm)*exp.collapsed_vector[i]*exp.collapsed_vector[j] << "\n";
			if(i==j){
				if (obs[i] >0 && exp[i] >0){
					cov[i][j] += 3.0 / (1.0/exp[i] + 2.0/obs[i]);
					// std::cout << i << " " << j << " " << 3.0 / (1.0/exp[i] + 2.0/obs[i]) << "\n";
					// poisson version
					// cov = (M-mu)**2 / (2*(mu - M + M*np.log(M/mu)))
					// cov[i][j] += pow(truespec_first.collapsed_vector[i]-thisGridCVSpec.collapsed_vector[i],2) / (2.0*(thisGridCVSpec.collapsed_vector[i]- truespec_first.collapsed_vector[i]+ fakeData[i]*log(fakeData[i]/truespec_second.collapsed_vector[i])));
				}
				else {
					cov[i][j] += exp[i]/2.0;
					// std::cout << i << " " << j << " " << exp[i]/2.0 << "\n";
				}
				// Old way
				// cov[i][j] += truespec_second[i];
			}
		}
	}
}

void dumpCoordFile(int dm2_grdpts, double dm2_lowbound, double dm2_hibound, int sin22th_grdpts, double sin22th_lowbound, bool isDecay){
	std::ofstream coordfile;
	double sin22th;
	double mnu;
  if (isDecay){
    coordfile.open("bins_decay.txt", std::ios_base::app);//std::ios_base::app

  }
  else{
    coordfile.open("bins.txt", std::ios_base::app);//std::ios_base::app
  }

	// Print binedges for easier plotting **
	for(int mi = 0; mi <= dm2_grdpts; mi++){
		mnu = pow(10.,((mi)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		coordfile << pow(mnu,2) << " ";
	}
	coordfile << std::endl;
	for(int si = 0; si <= sin22th_grdpts; si++){
		sin22th = pow(10.,(si/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
		coordfile << sin22th << " ";
	}
	coordfile << std::endl;
}

void getPrecompSpec(std::vector<sbn::SBNspec>& vecSpec, int dm2_grdpts, double dm2_hibound, double dm2_lowbound, std::string xml, double scaleFactor, bool islogname){
	std::cout <<"Precomputing the Mass Term Spectra\n";
	double mnu;
	for(int mi = 0; mi < dm2_grdpts; mi++){
		mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		std::stringstream stream;
    if (islogname){
      stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
    }
    else{
      stream << std::fixed << std::setprecision(4) << pow(mnu,2);
    }
		std::string infile = "DL_SINSQ_dm_"+stream.str()+".SBNspec.root";
		// if (scaleFactor != 1.0){
		// 	infile = "DL_SINSQ_dm_1x_"+stream.str()+".SBNspec.root";
		// }
		sbn::SBNspec inspec = sbn::SBNspec(infile,xml);
		inspec.RemoveMCError();
		inspec.ScaleAll(scaleFactor);
		// inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
		inspec.CollapseVector();
		vecSpec.push_back( inspec );
	}
}

void getPrecompSpecDecay(std::vector<std::vector<sbn::SBNspec>>& vecSpec, int dm2_grdpts, double dm2_hibound, double dm2_lowbound, int Um4sq_grdpts, double Um4sq_hibound, double Um4sq_lowbound, std::string xml, double scaleFactor){
	std::cout <<"Loading the Decay Spectra\n";
	double mnusq;
  double Um4;
	for(int mi = 0; mi < dm2_grdpts; mi++){
    std::vector<sbn::SBNspec> thismiVecSpec;
    for(int ui = 0;ui < Um4sq_grdpts;ui++){
      mnusq = pow(pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound)))),2);
      Um4 = sqrt(pow(10.,((ui+.5)/float(Um4sq_grdpts)*TMath::Log10(sqrt(Um4sq_hibound)/sqrt(Um4sq_lowbound)) + TMath::Log10(sqrt(Um4sq_lowbound)))));


      std::stringstream stream;
      stream << std::fixed << std::setprecision(4) << mnusq;

      std::stringstream stream2;
      stream2 << std::fixed << std::setprecision(4) << Um4;


      std::string infile = "DLDEC_um4_"+stream2.str()+"_dm_"+stream.str()+".SBNspec.root";
      // if (scaleFactor != 1.0){
        // 	infile = "DL_SINSQ_dm_1x_"+stream.str()+".SBNspec.root";
        // }
        sbn::SBNspec inspec = sbn::SBNspec(infile,xml);
        inspec.RemoveMCError();
        inspec.ScaleAll(scaleFactor);
        // inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
        inspec.CollapseVector();
        thismiVecSpec.push_back( inspec );
    }
    vecSpec.push_back( thismiVecSpec );
	}
}

void scaleDisappearedSpec(sbn::SBNspec& cvSpec, sbn::SBNspec massSpec, double disFactor, bool doPrint){
	std::vector<double> cvBefore(19,0);
	std::vector<double> oscBefore(19,0);
	std::vector<double> oscAfter(19,0);
	std::vector<double> cvAfter(19,0);
	if ((doPrint) && (disFactor > 0.9)){
		std::cout << "With Dis " << disFactor<< "\n";
		for (auto& h :cvSpec.hist){
			std::cout << h.GetName() << ",";
		}
		std::cout <<"\n";
		for (int i=0;i<19;i++){
			for (auto& h :cvSpec.hist){
				std::cout << h.GetBinContent(i+1) << ",";
			}
			std::cout <<"\n";
			// std::cout << cvSpec.collapsed_vector[i] << " ";
			cvBefore[i] = cvSpec.collapsed_vector[i];
		}
	}
	cvSpec.Scale("fullosc",0.0);
	if ((doPrint) && (disFactor > 0.9)){
		// for (auto& h : massSpec.hist){
		// 	std::cout << h.GetName() << "\n";
		// }
		for (int i=0;i<19;i++){
			oscBefore[i] = massSpec.collapsed_vector[i];
		}
	}
	massSpec.Scale("bnb",-1.0*disFactor);
	massSpec.Scale("nue",0.0);
	massSpec.Scale("ccpi0",-1.0*disFactor);
	massSpec.Scale("ncpi0",-1.0*disFactor);
	massSpec.Scale("ext",0.0);
	if ((doPrint) && (disFactor > 0.9)){
		for (int i=0;i<19;i++){
			oscAfter[i] = massSpec.collapsed_vector[i];
		}
	}
	// std::cout << cvSpec.xmlname << " " << massSpec.xmlname << "\n";
	cvSpec.Add(&massSpec);
	if ((doPrint) && (disFactor > 0.9)){
		for (int i=0;i<19;i++){
			cvAfter[i] = cvSpec.collapsed_vector[i];
		}
	}
	if ((doPrint) && (disFactor > 0.9)){
		std::cout <<"\n";
		for (int i=0;i<19;i++){
			std::cout << cvBefore[i] << "," << oscBefore[i] << "," << oscAfter[i] << "," << cvAfter[i] << "\n";
		}
	}
}



TMatrixD getZeroedMat(int nBins){
	TMatrixD covFracSys_second;
	covFracSys_second.ResizeTo(nBins,nBins);
	for (int i =0;i<nBins;i++){
		for (int j=0;j<nBins;j++){
			covFracSys_second[i][j] = 0.0;
		}
	}
	return covFracSys_second;
}

std::vector<std::vector<TMatrixD>> genCovMats(const sbn::SBNspec& cvSpec,
																							const std::vector<sbn::SBNspec>& a_sinsqSpec,
																							const TMatrixD& covFracSys,
																							const int dm2_grdpts,
																							const int sin22th_grdpts,
																							const double sin22th_lowbound,
																							const int nBins
																							){
	std::vector<std::vector<TMatrixD>> covMats(dm2_grdpts,std::vector<TMatrixD>(dm2_grdpts,TMatrixD(getZeroedMat(nBins))));
	for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
		for(int sin22thi_base = 0; sin22thi_base <sin22th_grdpts; sin22thi_base++){

			double sin22th_first = pow(10.,((sin22thi_base+0.5)/float(sin22th_grdpts)*TMath::Log10(1./sin22th_lowbound) + TMath::Log10(sin22th_lowbound)));
			sbn::SBNspec truespec_first = cvSpec;
			scaleDisappearedSpec(truespec_first, a_sinsqSpec[mi_base], sin22th_first);

			sbn::SBNchi TrueChi(truespec_first, covFracSys);
			TMatrixD covFracSys_base = getZeroedMat(nBins);
			TrueChi.FillCollapsedFractionalMatrix(&covFracSys_base);
			covMats[mi_base][sin22thi_base] = (covFracSys_base);
		}
	}
	return covMats;
}

void makeSpectraPlots(const sbn::SBNspec& cv, const sbn::SBNspec& app, const sbn::SBNspec& tot){
	gStyle->SetOptStat(0);
	TH1F cvSpec_h("cvSpec_h","cvSpec_h",19,250,1200);
	TH1F appSpec_h("appSpec_h","appSpec_h",19,250,1200);
	TH1F totSpec_h("totSpec_h","totSpec_h",19,250,1200);
	TH1F zeroSpec_h("zeroSpec_h","zeroSpec_h",19,250,1200);

	cvSpec_h.SetLineColor(kBlack);
	appSpec_h.SetLineColor(kRed);
	totSpec_h.SetLineColor(kBlue);
	zeroSpec_h.SetLineColor(kBlack);
	zeroSpec_h.SetLineStyle(kDashed);

	double maxVal = 0;
	double minVal = 999999;
	for (int i = 0; i<19; i++){
		if (cv.collapsed_vector[i] > maxVal) maxVal = cv.collapsed_vector[i];
		if (app.collapsed_vector[i] > maxVal) maxVal = app.collapsed_vector[i];
		if (tot.collapsed_vector[i] > maxVal) maxVal = tot.collapsed_vector[i];
		if (cv.collapsed_vector[i] < minVal) minVal = cv.collapsed_vector[i];
		if (app.collapsed_vector[i] < minVal) minVal = app.collapsed_vector[i];
		if (tot.collapsed_vector[i] < minVal) minVal = tot.collapsed_vector[i];
		cvSpec_h.SetBinContent(i+1,cv.collapsed_vector[i]);
		appSpec_h.SetBinContent(i+1,app.collapsed_vector[i]);
		totSpec_h.SetBinContent(i+1,tot.collapsed_vector[i]);
		zeroSpec_h.SetBinContent(i+1,0.0);
	}
	TCanvas c1("c","c",1200,1000);
	cvSpec_h.SetMaximum(maxVal*1.1);
	cvSpec_h.SetMinimum(minVal*1.1);
	cvSpec_h.Draw();
	appSpec_h.Draw("SAME");
	totSpec_h.Draw("SAME");
	zeroSpec_h.Draw("SAME");
	TLegend legend(0.7,0.70,0.90,0.9);
	legend.AddEntry(&cvSpec_h,"CV");
	legend.AddEntry(&appSpec_h,"Disap");
	legend.AddEntry(&totSpec_h,"Total");
	legend.Draw();
	c1.SaveAs("Spectras.png");
	assert (1==2);
}

void makeSpectraPlotsComp(const std::vector<double>& cv, const std::vector<double>& disap, int sin_base, int mi_base, double sin_val, double mi_val, int nBins, std::string folder){
	gStyle->SetOptStat(0);
	std::string titleStr = "#nu_{mu} Dis. Sin^2 Idx: "+std::to_string(sin_base)+", dM Idx: " +std::to_string(mi_base);
	TH1F cvSpec_h(titleStr.c_str(),titleStr.c_str(),nBins,250,1200);
	TH1F disSpec_h("tmp","tmp",nBins,250,1200);
  // TH1F cvSpec_h(titleStr.c_str(),titleStr.c_str(),nBins,0,600000);
  // TH1F disSpec_h("tmp","tmp",nBins,0,600000);
	cvSpec_h.SetLineColor(kBlack);
	disSpec_h.SetLineColor(kBlue);
  cvSpec_h.SetLineWidth(5);
  disSpec_h.SetLineWidth(5);
	// disSpec_h.SetLineStyle(kDashed);
	double maxVal = 0;
	double minVal = 0;
	for (int i = 0; i<nBins; i++){
		if (cv[i] > maxVal) maxVal = cv[i];
		// if (disap[i] > maxVal) maxVal = disap[i];
		if (cv[i] < minVal) minVal = cv[i];
		// if (disap[i] < minVal) minVal = disap[i];
		cvSpec_h.SetBinContent(i+1,cv[i]);
		disSpec_h.SetBinContent(i+1,disap[i]);
	}
	TCanvas c1("c","c",1200,1000);
  cvSpec_h.SetXTitle("Neutrino Energy");
  // cvSpec_h.SetXTitle("Q2 Reco");

  cvSpec_h.SetYTitle("Events");
	cvSpec_h.SetMaximum(maxVal*1.1);
  disSpec_h.SetMaximum(maxVal*1.1);
	cvSpec_h.Draw();
	disSpec_h.Draw("SAME");
	TLegend legend(0.45,0.75,0.90,0.9);
	std::string cv_str  = "Null Model Spectra";
	std::string dis_str = "Dis. Sin^2: "+to_string_with_precision(std::round(sin_val*1000.0)/1000) + " dM^2: " +to_string_with_precision(std::round(mi_val*1000.0)/1000);
  legend.AddEntry(&cvSpec_h, cv_str.c_str());
	legend.AddEntry(&disSpec_h,dis_str.c_str());
	legend.Draw();
	std::string savestr = folder + "/spec_m"+ZeroPadNumber(mi_base,2) + "_s"+ZeroPadNumber(sin_base,2)+".png";
	c1.SaveAs(savestr.c_str());
}


void makeSpectraPlotsComp_Ratio(const std::vector<double>& cv, const std::vector<double>& disap, int sin_base, int mi_base, double sin_val, double mi_val, int nBins, std::string folder){
	gStyle->SetOptStat(0);
	// std::string titleStr = "#nu_{#mu} Dis. Sin^2 Idx: "+std::to_string(sin_base)+", dM Idx: " +std::to_string(mi_base);
  std::string titleStr = "Muon Neutrino Disappearance";
  TH1F cvSpec_h(titleStr.c_str(),titleStr.c_str(),nBins,250,1200);
  TH1F disSpec_h("tmp","tmp",nBins,250,1200);
  TH1F fakedata_h("tmp2","tmp2",nBins,250,1200);
  // TH1F cvSpec_h(titleStr.c_str(),titleStr.c_str(),nBins,0,600000);
  // TH1F disSpec_h("tmp","tmp",nBins,0,600000);
  // TH1F fakedata_h("tmp2","tmp2",nBins,0,600000);
	cvSpec_h.SetLineColor(kBlack);
	disSpec_h.SetLineColor(kBlue);
  fakedata_h.SetLineColor(kRed);
  cvSpec_h.SetLineWidth(5);
  disSpec_h.SetLineWidth(5);
  fakedata_h.SetLineWidth(5);
  // std::vector<double> f_v = {35,251,402,481,433,401,290,227,179,154,132,97,75,55};
	// disSpec_h.SetLineStyle(kDashed);
	double maxVal = 0;
	double minVal = 0;
	for (int i = 0; i<nBins; i++){
		if (cv[i] > maxVal) maxVal = cv[i];
		// if (disap[i] > maxVal) maxVal = disap[i];
		if (cv[i] < minVal) minVal = cv[i];
		// if (disap[i] < minVal) minVal = disap[i];
		cvSpec_h.SetBinContent(i+1,cv[i]);
		disSpec_h.SetBinContent(i+1,disap[i]);
    // fakedata_h.SetBinContent(i+1,f_v[i]);
	}
	// TCanvas c1("c","c",1200,1000);
  // TPad  pad("upper_pad","",0,0.26,1.,1.);
  // TPad rpad("ratio_pad","",0,0,1.,0.24);
  // pad.SetFrameLineWidth(2);
  // pad.SetTopMargin(.1);
  // pad.SetBottomMargin(0.0);
  // pad.SetLeftMargin(0.1);
  // rpad.SetFrameLineWidth(2);
  // rpad.SetTopMargin(0.05);
  // rpad.SetBottomMargin(0.4);
  // rpad.SetLeftMargin(0.1);

  TCanvas c1("c","c",1200,1000);
  double botfrac = 0.23;
  double topfrac = botfrac + 0.01;
  double rightmargin = 0.12;
  double leftmargin = 0.12;

  TPad  pad("upper_pad","",0,topfrac,1.,1.);
  TPad rpad("ratio_pad","",0,0,1.,botfrac);
  pad.SetFrameLineWidth(2);
  pad.SetTopMargin(0.88*1/topfrac);
  pad.SetBottomMargin(0.0);
  pad.SetLeftMargin(leftmargin);
  pad.SetRightMargin(rightmargin);
  rpad.SetFrameLineWidth(2);
  rpad.SetTopMargin(0.05);
  rpad.SetBottomMargin(0.12*1/botfrac);
  rpad.SetLeftMargin(leftmargin);
  rpad.SetRightMargin(rightmargin);

  // gStyle.SetPadTickX(1)
  // gStyle.SetPadTickY(1)
  // gStyle.SetLineWidth(2)
  pad.Draw();
  pad.cd();
  // TH1F frame = pad.DrawFrame(xmin,ymin,xmax,ymax)
  c1.cd();
  pad.cd();


  cvSpec_h.SetXTitle("Neutrino Energy");
  // cvSpec_h.SetXTitle("Q2 Reco");

  cvSpec_h.SetYTitle("Events");
	cvSpec_h.SetMaximum(maxVal*1.1);
  disSpec_h.SetMaximum(maxVal*1.1);
	cvSpec_h.Draw();
	disSpec_h.Draw("SAME");
  // fakedata_h.Draw("SAME");
	TLegend legend(0.45,0.75,1-rightmargin,0.9);
	std::string cv_str  = "Null Model Spectra";
	// std::string dis_str = "Dis. Sin^2: "+to_string_with_precision(std::round(sin_val*1000.0)/1000) + " dM^2: " +to_string_with_precision(std::round(mi_val*1000.0)/1000);
  std::string dis_str = "Dis. Sin^{2}2#theta: "+to_string_with_precision(std::round(sin_val*1000.0)/1000) + " #Deltam^{2}: " +to_string_with_precision(std::round(mi_val*1000.0)/1000);
  // std::string dis_str = "Dis. U^{2}_{#mu4}: "+to_string_with_precision(std::round(sin_val*1000.0)/1000) + " #Deltam^{2}: " +to_string_with_precision(std::round(mi_val*1000.0)/1000);

  legend.AddEntry(&cvSpec_h, cv_str.c_str());
	legend.AddEntry(&disSpec_h,dis_str.c_str());
  // legend.AddEntry(&fakedata_h,"Fake Dataset 1");
	legend.Draw();


  c1.cd();
  rpad.SetBorderSize(1);
  rpad.Draw();
  rpad.cd();

  // TH1F rframe = (*rpad.DrawFrame(0,0.0,600000,1.2));
  TH1F rframe = (*rpad.DrawFrame(250,0.0,1200,1.2));

  rframe.GetYaxis()->SetTitle("Osc. / Null");
  rframe.GetXaxis()->SetLabelFont(43);
  rframe.GetXaxis()->SetLabelSize(25);
  rframe.GetYaxis()->SetLabelFont(43);
  rframe.GetYaxis()->SetLabelSize(25);
  rframe.GetYaxis()->SetNdivisions(7);
  rframe.GetXaxis()->SetTitleFont(43);
  rframe.GetXaxis()->SetTitleSize(25);
  rframe.GetYaxis()->SetTitleFont(43);
  rframe.GetYaxis()->SetTitleSize(25);
  rframe.GetYaxis()->SetTitleOffset(1.5);
  // rframe.GetXaxis()->SetTitle("Q2 Reco");
  rframe.GetXaxis()->SetTitle("Neutrino Energy");
  rframe.GetXaxis()->SetTitleOffset(5);
  rframe.Draw();
  // TH1F cvSpec_ratio("tmp3","tmp3",nBins,0,60000);
  // TH1F disSpec_ratio("tmp4","tmp4",nBins,0,60000);
  // TH1F fakedata_ratio("tmp5","tmp5",nBins,0,60000);
  TH1F cvSpec_ratio("tmp3","tmp3",nBins,250,1200);
  TH1F disSpec_ratio("tmp4","tmp4",nBins,250,1200);
  TH1F fakedata_ratio("tmp5","tmp5",nBins,250,1200);
  cvSpec_ratio.SetLineColor(kBlack);
  disSpec_ratio.SetLineColor(kBlue);
  fakedata_ratio.SetLineColor(kRed);
  cvSpec_ratio.SetLineWidth(5);
  disSpec_ratio.SetLineWidth(5);
  fakedata_ratio.SetLineWidth(5);
  for (int i=0;i<nBins;i++){
    cvSpec_ratio.SetBinContent(i+1,1.0);
    disSpec_ratio.SetBinContent(i+1,disSpec_h.GetBinContent(i+1)*1.0/cvSpec_h.GetBinContent(i+1));
    fakedata_ratio.SetBinContent(i+1,fakedata_h.GetBinContent(i+1)*1.0/cvSpec_h.GetBinContent(i+1));
  }
  cvSpec_ratio.Draw("SAME");
  disSpec_ratio.Draw("SAME");
  // fakedata_ratio.Draw("SAME");

	std::string savestr = folder + "/spec_m"+ZeroPadNumber(mi_base,2) + "_s"+ZeroPadNumber(sin_base,2)+".png";
	c1.SaveAs(savestr.c_str());
}


void makeSpectraPlotsComp_Data(const std::vector<double>& cv, const std::vector<double>& disap, const std::vector<double>& best, int set, const std::vector<double>& unc){
	gStyle->SetOptStat(0);
  int nBins = 19;
	// std::string titleStr = "#nu_{mu} Dis. Fake Dataset: "+std::to_string(set);
  std::string titleStr = "Muon Neutrino Disappearance Fake Dataset: "+std::to_string(set);
  if (set == 3){
    titleStr = "DL Muon Neutrino Disappearance, Data: 6.67e20 POT";
  }
  TH1F label_h(titleStr.c_str(),titleStr.c_str(),19,250,1200);
	TH1F cvSpec_h("tmp0","tmp0",19,250,1200);
	TH1F disSpec_h("tmp","tmp",19,250,1200);
  TH1F best_h("tmp2","tmp2",19,250,1200);

	cvSpec_h.SetLineColor(kBlue);
	disSpec_h.SetLineColor(kBlack);
  best_h.SetLineColor(kRed);
  cvSpec_h.SetLineWidth(5);
  disSpec_h.SetLineWidth(5);
  best_h.SetLineWidth(5);
  disSpec_h.SetMarkerStyle(20);


  TGraphErrors errors_cv;
  errors_cv.SetLineColor(kBlue-9);
  errors_cv.SetLineWidth(2.5);
  errors_cv.SetFillStyle(3353);
  gStyle->SetHatchesLineWidth(2.5);

  errors_cv.SetFillColor(kBlue-9);

	// disSpec_h.SetLineStyle(kDashed);
	double maxVal = 0;
	double minVal = 0;
	for (int i = 0; i<19; i++){
		if (cv[i] > maxVal) maxVal = cv[i];
		// if (disap[i] > maxVal) maxVal = disap[i];
		if (cv[i] < minVal) minVal = cv[i];
		// if (disap[i] < minVal) minVal = disap[i];
		cvSpec_h.SetBinContent(i+1,cv[i]);
		disSpec_h.SetBinContent(i+1,disap[i]);
    best_h.SetBinContent(i+1,best[i]);

    disSpec_h.SetBinError(i+1,pow(disap[i],0.5));
    errors_cv.SetPoint(i,275+i*50,cv[i]);
    errors_cv.SetPointError(i,25,unc[i]);

    // cvSpec_h.SetBinError(i+1, unc[i]);
	}
  TCanvas c1("c","c",1200,1000);
  double botfrac = 0.23;
  double topfrac = botfrac + 0.01;
  double rightmargin = 0.12;
  double leftmargin = 0.12;

  TPad  pad("upper_pad","",0,topfrac,1.,1.);
  TPad rpad("ratio_pad","",0,0,1.,botfrac);
  pad.SetFrameLineWidth(2);
  pad.SetTopMargin(0.88*1/topfrac);
  pad.SetBottomMargin(0.0);
  pad.SetLeftMargin(leftmargin);
  pad.SetRightMargin(rightmargin);
  rpad.SetFrameLineWidth(2);
  rpad.SetTopMargin(0.05);
  rpad.SetBottomMargin(0.12*1/botfrac);
  rpad.SetLeftMargin(leftmargin);
  rpad.SetRightMargin(rightmargin);

  pad.Draw();
  pad.cd();
  // TH1F frame = pad.DrawFrame(xmin,ymin,xmax,ymax)
  c1.cd();
  pad.cd();

  label_h.SetXTitle("Reconstructed Neutrino Energy (MeV)");
  label_h.SetYTitle("Events");
	label_h.SetMaximum(maxVal*1.2);
  disSpec_h.SetMaximum(maxVal*1.2);
  label_h.Draw();
  errors_cv.Draw("SAME 2");
	cvSpec_h.Draw("SAME");

	disSpec_h.Draw("SAMEP");
  best_h.Draw("SAME");

  TLegend legend(0.45,0.75,1-rightmargin,0.9);
	std::string cv_str  = "Null Model Spectra";
	std::string dis_str = "Fake Dataset "+std::to_string(set);
  std::string best_sr = "Best Fit Disappearance Spectra";
  if (set ==3){
    cv_str  = "Null Oscillation Model";
    dis_str = "6.67e20 Data";
    best_sr = "Best Fit 3+1 Disappearance Model";
  }
  legend.AddEntry(&cvSpec_h, cv_str.c_str());
	legend.AddEntry(&disSpec_h,dis_str.c_str());
  legend.AddEntry(&best_h,best_sr.c_str());
	legend.Draw();

  c1.cd();
  rpad.SetBorderSize(1);
  rpad.Draw("SAME");
  rpad.cd();
  double minrat =0.0;
  double maxrat =1.2;
  if (set == 3){
    minrat = 0.6;
    maxrat = 1.4;
  }
  TH1F rframe = (*rpad.DrawFrame(250,minrat,1200,maxrat));

  rframe.GetYaxis()->SetTitle("Ratio to Null");
  rframe.GetXaxis()->SetLabelFont(43);
  rframe.GetXaxis()->SetLabelSize(25);
  rframe.GetYaxis()->SetLabelFont(43);
  rframe.GetYaxis()->SetLabelSize(25);
  rframe.GetYaxis()->SetNdivisions(7);
  if (set == 3){
    rframe.GetYaxis()->SetNdivisions(5);
  }
  rframe.GetXaxis()->SetTitleFont(43);
  rframe.GetXaxis()->SetTitleSize(25);
  rframe.GetYaxis()->SetTitleFont(43);
  rframe.GetYaxis()->SetTitleSize(25);
  rframe.GetYaxis()->SetTitleOffset(1.5);
  // rframe.GetXaxis()->SetTitle("Q2 Reco");
  rframe.GetXaxis()->SetTitle("Reconstructed Neutrino Energy (MeV)");
  rframe.GetXaxis()->SetTitleOffset(5);
  rframe.Draw();
  // TH1F cvSpec_ratio("tmp3","tmp3",nBins,0,60000);
  // TH1F disSpec_ratio("tmp4","tmp4",nBins,0,60000);
  // TH1F fakedata_ratio("tmp5","tmp5",nBins,0,60000);
  TH1F label_ratio("tmp1","tmp1",nBins,250,1200);
  TH1F cvSpec_ratio("tmp3","tmp3",nBins,250,1200);
  TH1F best_ratio("tmp4","tmp4",nBins,250,1200);
  TH1F disap_ratio("tmp5","tmp5",nBins,250,1200);
  cvSpec_ratio.SetLineColor(kBlue);
  best_ratio.SetLineColor(kRed);
  disap_ratio.SetLineColor(kBlack);
  cvSpec_ratio.SetLineWidth(5);
  best_ratio.SetLineWidth(5);
  disap_ratio.SetLineWidth(5);
  disap_ratio.SetMarkerStyle(20);


  TGraphErrors errors_rat_cv;
  errors_rat_cv.SetLineColor(kBlue-9);
  errors_rat_cv.SetLineWidth(2.5);
  errors_rat_cv.SetFillStyle(3353);
  errors_rat_cv.SetFillColor(kBlue-9);

  for (int i=0;i<nBins;i++){
    cvSpec_ratio.SetBinContent(i+1,1.0);
    best_ratio.SetBinContent(i+1,best_h.GetBinContent(i+1)*1.0/cvSpec_h.GetBinContent(i+1));
    disap_ratio.SetBinContent(i+1,disSpec_h.GetBinContent(i+1)*1.0/cvSpec_h.GetBinContent(i+1));

    errors_rat_cv.SetPoint(i,275+i*50,1.0);
    errors_rat_cv.SetPointError(i,25,unc[i]/cvSpec_h.GetBinContent(i+1));
    disap_ratio.SetBinError(i+1,pow(disap[i],0.5)/disap[i]);

    // cvSpec_ratio.SetBinError(i+1, unc[i]/cvSpec_h.GetBinContent(i+1));
  }
  label_ratio.Draw("SAME");
  errors_rat_cv.Draw("SAME 2");
  cvSpec_ratio.Draw("SAME");
  best_ratio.Draw("SAME");
  disap_ratio.Draw("SAMEP");

	std::string savestr = "data_spectras/spec_d"+std::to_string(set)+".png";
	c1.SaveAs(savestr.c_str());
}


sbn::SBNspec getAppSpec_exact(float mnu, bool doGen){
	std::string tag = "DL";
	std::string xml = "/cluster/tufts/wongjiradlab/jmills09/whipping_star/xml/numudisappearance.xml";
	if (doGen == true){
		// if(false){
		// 		sbn::SBNcovariance _covar(xml);
		// 	 _covar.FormCovarianceMatrix(tag);
		// 	 _covar.PrintMatricies(tag);
		// 	 _covar.frac_covariance.Print();
		// }
		// start with a null model
		sbn::NeutrinoModel nullModel(0, 0, 0);
		sbn::SBNgenerate * bkgo = new sbn::SBNgenerate(xml,nullModel);
		sbn::SBNspec bkg = bkgo->spec_central_value;
		bkg.RemoveMCError();
		//SBNspec bkg("DL.SBNspec.root",xml);
		// set full osc to 0 for bkg spectrum
		// bkg.Scale("fullosc",0.0);
		bkg.WriteOut(tag+"_Bkg");
		//Model: mnu, ue4, um4
		sbn::NeutrinoModel testModel(mnu, sqrt(.5), sqrt(.5)); //NuMu Disap Maximized
		// on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
		sbn::SBNgenerate * gen = new sbn::SBNgenerate(xml,testModel);
		// Write them to files
		gen->WritePrecomputedOscSpecs(tag);
	}

	std::stringstream stream;
	stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
	std::string infile = "DL_SINSQ_dm_"+stream.str()+".SBNspec.root";
	sbn::SBNspec inspec = sbn::SBNspec(infile,xml);
	inspec.RemoveMCError();
	// inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
	inspec.CollapseVector();
	return inspec;
}
