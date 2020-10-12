#include "SBNsinglephoton.h"

using namespace sbn;


SBNsinglephoton::SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid):SBNsinglephoton(xmlname, intag, ingrid,ingrid, false){}

SBNsinglephoton::SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid, NGrid in_polygrid, bool has_polygrid): SBNconfig(xmlname), tag(intag), m_grid(ingrid), m_bool_poly_grid(has_polygrid){

    if(is_verbose) std::cout << "SBNsinglephoton::SBNsinglephoton\t|| Setup grid" << std::endl;
    m_vec_grid = m_grid.GetGrid();
    m_num_total_gridpoints = m_grid.f_num_total_points;

    m_poly_total_gridpoints = 1;
    if(m_bool_poly_grid){
	 m_poly_grid=in_polygrid;
	 m_vec_poly_grid = m_poly_grid.GetGrid();
         m_poly_total_gridpoints = m_poly_grid.f_num_total_points;
    }

    m_total_gridpoints = m_num_total_gridpoints*m_poly_total_gridpoints;
    m_max_number_iterations = 20;
    m_chi_min_convergance_tolerance = 0.001;
    m_bool_modify_cv = false;
    m_bool_cv_spectrum_generated = false;
    m_bool_cv_spectrum_loaded = false;
    m_bool_data_spectrum_loaded = false;
    m_interpolation_number=100;

    m_full_fractional_covariance_matrix = NULL;
    m_full_but_genie_fractional_covariance_matrix = NULL;
    m_genie_fractional_covariance_matrix = NULL;

    m_file_open=false;
}

int SBNsinglephoton::SetPolyGrid(NGrid in_polygrid){
    if(is_verbose) std::cout << "SBNsinglephoton::SetPolyGrid\t|| Setup polynomial grid " << std::endl;
    m_bool_poly_grid = true;
    m_poly_grid=in_polygrid;
    m_vec_poly_grid = m_poly_grid.GetGrid();
    m_poly_total_gridpoints = m_poly_grid.f_num_total_points;

    m_total_gridpoints = m_num_total_gridpoints*m_poly_total_gridpoints;
 
    return 0;
}

int SBNsinglephoton::OpenFiles(){

    num_files = montecarlo_file.size();
    montecarlo_additional_weight.resize(num_files,1.0);
    montecarlo_additional_weight_formulas.resize(num_files);


    if(is_verbose) std::cout<< "SBNsinglephoton::OpenFiles\t|| Open the files"<< std::endl;
    for(auto &fn: montecarlo_file){
        files.push_back(new TFile(fn.c_str()));
        if(files.back()->IsZombie() || !files.back()->IsOpen()){
            std::cout<<"SBNsinglephoton::OpenFiles\t|| ERROR! Failed to open the file "<<fn<<std::endl;
            exit(EXIT_FAILURE);
        }
    }


    for(int i=0; i<montecarlo_name.size(); i++){
        std::cout<<"Getting TTree "<<montecarlo_name[i]<<" from file "<<montecarlo_file[i]<<std::endl;
        trees.push_back((TTree*)files.at(i)->Get(montecarlo_name.at(i).c_str()) );
        std::cout<<"--TTree has "<<trees.back()->GetEntries()<<" entries. "<<std::endl;
    }

   if(is_verbose) std::cout << "SBNsinglephoton::OpenFiles\t||Setup any friendtrees" << std::endl;
   for(int i=0; i<montecarlo_file.size(); i++){
        const auto& fn = montecarlo_file.at(i);
        auto montecarlo_file_friend_treename_iter = montecarlo_file_friend_treename_map.find(fn);
        if (montecarlo_file_friend_treename_iter != montecarlo_file_friend_treename_map.end()) {
            std::cout<<" Detected friend trees "<<std::endl;

            auto montecarlo_file_friend_iter = montecarlo_file_friend_map.find(fn);
            if (montecarlo_file_friend_iter == montecarlo_file_friend_map.end()) {
                std::stringstream ss;
                ss << "Looked for filename=" << fn << " in fnmontecarlo_file_friend_iter, but could not be found... bad config?" << std::endl;
                throw std::runtime_error(ss.str());
            }

            for(int k=0; k < (*montecarlo_file_friend_iter).second.size(); k++){
                std::string treefriendname = (*montecarlo_file_friend_treename_iter).second.at(k);
                std::string treefriendfile = (*montecarlo_file_friend_iter).second.at(k);

                std::cout <<" Adding a friend tree:  " <<treefriendname<<" from file: "<< treefriendfile <<std::endl;

                trees[i]->AddFriend(treefriendname.c_str(),treefriendfile.c_str());
            }
        }

    }

    for(auto &t: trees){
        nentries.push_back(t->GetEntries());
    }


    for(int i=0; i< num_files; i++){

        double pot_scale = 1.0;
        if(montecarlo_pot[i]!=-1){
            pot_scale = this->plot_pot/montecarlo_pot[i];
        }

        montecarlo_scale[i] = montecarlo_scale[i]*pot_scale;

        std::cout << " TFile::Open() file=" << files[i]->GetName() << " @" << files[i] << std::endl;
        std::cout << " Has POT " <<montecarlo_pot[i] <<" and "<<nentries[i] <<" entries "<<std::endl;

        for(int k=0; k<branch_variables.at(i).size(); k++){
            const auto branch_variable = branch_variables.at(i).at(k);
            std::cout<<"Setting Branch: "<<branch_variable->name<<std::endl;
             branch_variable->branch_formula =  new TTreeFormula(("branch_form"+std::to_string(i)).c_str(), branch_variable->name.c_str(), trees[i]);

            if(branch_variable->GetOscillate()){
                std::cout<<"Setting true branch variables"<<std::endl;
		//if true_param_name is given as a branch name
                //trees.at(i)->SetBranchAddress( branch_variable->true_param_name.c_str(), branch_variable->GetTrueValue() );

		//if true_param_name is given in the form of a formula
		branch_variable->branch_true_value_formula = new TTreeFormula(("branch_true_value_form"+std::to_string(i)).c_str(),branch_variable->true_param_name.c_str(), trees[i]);
            }
        }

        if(montecarlo_additional_weight_bool[i]){
            std::cout<<"Setting Additional weight of : "<< montecarlo_additional_weight_names[i].c_str()<<std::endl;
            montecarlo_additional_weight_formulas[i] =  new TTreeFormula(("a_w"+std::to_string(i)).c_str(),montecarlo_additional_weight_names[i].c_str(),trees[i]);
        }


    }

    m_file_open= true;
    if(is_verbose) std::cout << "SBNsinglephoton::OpenFiles\t|| Finish opening files and setting up TTrees " << std::endl; 

    return 0;
}


int SBNsinglephoton::CloseFiles(){
    if(is_verbose) std::cout<< "SBNsinglephoton::CloseFiles\t|| Closing TFiles..."<< std::endl;
    for(auto f: files){
            std::cout <<" TFile::Close() file=" << f->GetName() << " @" << f << std::endl;
            f->Close();
    }
    m_file_open=false;
    return 0;
}


int SBNsinglephoton::ScaleSpectrum(SBNspec* inspec, double flat_factor, std::vector<double>& param){
//    inspec = new SBNspec(xmlname,i,false);

    SBNspec temp = *inspec;
    //if(is_verbose) std::cout<<"SBNsinglephoton::ScaleSpectrum\t|| -----------------------------------------------\n";

    //std::cout << num_files <<" " <<  __LINE__ << std::endl;
    for(int j=0;j<num_files;j++){


        for(int i=0; i< std::min(  montecarlo_maxevents.at(j)  ,nentries.at(j)); i++){
            trees.at(j)->GetEntry(i);

            //if(i%1000==0) std::cout<<"SBNsinglephoton::ScaleSpectrum\t|| On event: "<<i<<" of "<<nentries[j]<<" from File: "<<montecarlo_file[j]<<std::endl;

            double global_weight = 1.0;
            if( montecarlo_additional_weight_bool[j]){
                    montecarlo_additional_weight_formulas[j]->GetNdata();
                    global_weight = montecarlo_additional_weight_formulas[j]->EvalInstance();
            };//this will be 1.0 unless specified
            global_weight = global_weight*montecarlo_scale[j];



            if(std::isinf(global_weight) || global_weight != global_weight){
                std::cout<<"SBNsinglephoton::ScaleSpectrum\t|| ERROR  error @ "<<i<<" in File "<<montecarlo_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
                exit(EXIT_FAILURE);
            }


            for(int t=0; t<branch_variables[j].size();t++){
                    const auto branch_variable = branch_variables[j][t];
                    int ih = inspec->map_hist.at(branch_variable->associated_hist);

		    //grab reco info and determine which reco bin it belongs to
                    branch_variable->GetFormula()->GetNdata();
                    double reco_var = branch_variable->GetFormula()->EvalInstance();
                    int reco_bin = inspec->GetGlobalBinNumber(reco_var,ih);

                    if(branch_variables[j][t]->GetOscillate()){
			//truth info
  			double true_var = *(static_cast<double*>(branch_variables[j][t]->GetTrueValue()));

                        double prescale_factor = this->ScaleFactor(true_var, flat_factor, param);

                        inspec->hist[ih].Fill(reco_var, global_weight*prescale_factor);
 		    }
            }
        } //end of entry loop
    } // end of nfile loop

    inspec->CalcFullVector();
    return 0;
}

int SBNsinglephoton::PreScaleSpectrum(std::string xmlname, double flat_factor, std::vector<double>& param){
    SBNspec tm(xmlname,-1,false);
    SBNspec temp_cv_spectrum = tm;
    SBNspec spec_prescale  = tm;   //prescaled spectrum

    std::cout<<"SBNsinglephoton::PreScaleSpectrum\t|| -----------------------------------------------\n";
    std::cout<<"SBNsinglephoton::PreScaleSpectrum\t|| -----------------------------------------------\n";

    for(int j=0;j<num_files;j++){


        for(int i=0; i< std::min(  montecarlo_maxevents.at(j)  ,nentries.at(j)); i++){
            trees.at(j)->GetEntry(i);

            if(i%100==0) std::cout<<"SBNsinglephoton::PreScaleSpectrum\t|| On event: "<<i<<" of "<<nentries[j]<<" from File: "<<montecarlo_file[j]<<std::endl;

            double global_weight = 1.0;
            if( montecarlo_additional_weight_bool[j]){
                    montecarlo_additional_weight_formulas[j]->GetNdata();
                    global_weight = montecarlo_additional_weight_formulas[j]->EvalInstance();
            };//this will be 1.0 unless specified
            global_weight = global_weight*montecarlo_scale[j];



            if(std::isinf(global_weight) || global_weight != global_weight){
                std::cout<<"SBNsinglephoton::PreScaleSpectrum\t|| ERROR  error @ "<<i<<" in File "<<montecarlo_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
                exit(EXIT_FAILURE);
            }



            for(int t=0; t<branch_variables[j].size();t++){
                    const auto branch_variable = branch_variables[j][t];
                    int ih = temp_cv_spectrum.map_hist.at(branch_variable->associated_hist);

		    //grab reco info and determine which reco bin it belongs to
                    branch_variable->GetFormula()->GetNdata();
                    double reco_var = branch_variable->GetFormula()->EvalInstance();
                    int reco_bin = temp_cv_spectrum.GetGlobalBinNumber(reco_var,ih);

                    if(branch_variables[j][t]->GetOscillate()){
			//truth info
  			double true_var = *(static_cast<double*>(branch_variables[j][t]->GetTrueValue()));

                        double prescale_factor = this->ScaleFactor(true_var, flat_factor, param);

                        spec_prescale.hist[ih].Fill(reco_var, global_weight*prescale_factor);
         		if(!m_bool_cv_spectrum_generated)  temp_cv_spectrum.hist[ih].Fill(reco_var,global_weight);
 		    }else{
                        if(!m_bool_cv_spectrum_generated) temp_cv_spectrum.hist[ih].Fill(reco_var,global_weight);
  		    }
            }
        } //end of entry loop
    } // end of nfile loop

 
    if(is_verbose) std::cout<< "SBNsinglephoton::PreScaleSpectrum\t||Write out spectra" << std::endl;
    std::ostringstream prescale_tag;
    if(flat_factor != 0) prescale_tag << "Flat_" << std::fixed<< std::setprecision(3) << flat_factor << "_PreScaled";
    else prescale_tag << "PreScaled";
    //generate tag for prescaled spectra
    for(int i=0; i< param.size(); i++){
	prescale_tag << "_" << std::fixed<< std::setprecision(3) << param[i]; 
    }   
    if(param.size() != 0) spec_prescale.WriteOut(tag+"_"+prescale_tag.str());
    //if(param.size() != 0) spec_prescale.WriteOut(tag+"_PreScaled"+prescale_tag.str());

    if(!m_bool_cv_spectrum_generated){
	   m_bool_cv_spectrum_generated = true;
	   temp_cv_spectrum.WriteOut(tag+"_CV");
    }
    return 0;
}


int SBNsinglephoton::PreScaleSpectrum(std::string xmlname, std::vector<double>& param){
	return this->PreScaleSpectrum(xmlname, 0.0, param);
}


int SBNsinglephoton::GeneratePreScaledSpectra(){
    if(!m_file_open) this->OpenFiles();
    if(!m_bool_poly_grid){
	std::cout << "SBNsinglephoton::GeneratePreScaledSpectra\t|| No grid for polynomial scaling present------------------\n"<< std::endl;
	std::cout << "SBNsinglephoton::GeneratePreScaledSpectra\t|| Gonna generate only CV, then exiting GeneratePreScaledSpectra ..." << std::endl;
	std::vector<double> empty{};
	this->PreScaleSpectrum(xmlname, empty);
    }else{
	
    	//m_vec_poly_grid = m_poly_grid.GetGrid();
    	//m_poly_total_gridpoints = m_poly_grid.f_num_total_points;

    	for(int i=0; i< m_poly_total_gridpoints; i++){
		std::vector<double> point = m_vec_poly_grid[i];
		this->PreScaleSpectrum(xmlname, point);	
    	}
    }

    this->CloseFiles();

    return 0;
}


int SBNsinglephoton::LoadSpectraOnTheFly(){
	if(is_verbose) std::cout << "SBNsinglephoton::LoadSpectraOnTheFly\t||\tCalculate and Load scaled spectra on the fly" << std::endl;

	if(!m_file_open) this->OpenFiles();
	int flat_grid_index=-1;
	int flat_grid_npoint=0;
	int period=0;
	for(int i=0; i<m_grid.f_num_dimensions; i++){
	    period *= m_grid.f_dimensions[i].f_N;
	    if(m_grid.f_dimensions[i].f_name == m_poly_grid.f_dimensions[0].f_name){
		 flat_grid_index=i;
		 flat_grid_npoint=m_grid.f_dimensions[i].f_N;
		 period =1;
	    }
	}

	//calculate the flat+poly part
	int vec_index=0;
	if(is_verbose) std::cout << "SBNsinglephoton::LoadSpectraOnTheFly\t||\tGenerating flat+poly part of the spectra"<< std::endl;
	SBNspec tm(xmlname, -1, false);
	std::vector<SBNspec> m_prescale_spec(m_poly_total_gridpoints*flat_grid_npoint, tm);	
	for(int i=0; i<m_poly_total_gridpoints; i++){
	    std::vector<double> ipoint = m_vec_poly_grid[i];
	    std::vector<double> fgrid_point = m_grid.f_dimensions[flat_grid_index].f_points;
	    for(int j=0; j< flat_grid_npoint; j++){
		if(is_verbose && (vec_index%50==0)) std::cout << "On Point " << vec_index << std::endl;
		this->ScaleSpectrum(&m_prescale_spec[vec_index], fgrid_point[j], ipoint);
		vec_index++;
	    }
	}

	vec_index=0;
	//get the total scaled spectra
	if(is_verbose) std::cout << "SBNsinglephoton::LoadSpectraOnTheFly\t||\tGrab the whole scaled spectra " <<std::endl;
	m_scaled_spec_grid.resize(m_total_gridpoints);
	for(int i=0; i<m_poly_total_gridpoints; i++){
	    for(int j=0; j<m_num_total_gridpoints; j++){

		if(vec_index% 1000 == 0) std::cout << "On Point "<< vec_index << std::endl;
		
		std::vector<double> jpoint = m_vec_grid[j];

		m_scaled_spec_grid[vec_index]=*m_cv_spectrum;
		for(int k=0; k<jpoint.size(); k++){
		      if(k == flat_grid_index) m_scaled_spec_grid[vec_index].Scale(m_grid.f_dimensions[k].f_name,0.0);
		      m_scaled_spec_grid[vec_index].Scale(m_grid.f_dimensions[k].f_name, jpoint[k]);
		}

		int flat_index = (j%(period*flat_grid_npoint))/period; 
		m_scaled_spec_grid[vec_index].Add(&m_prescale_spec[i*flat_grid_npoint+flat_index]);
		if(tag == "NCpi0") m_scaled_spec_grid[vec_index].Scale("NCDeltaRadOverlayLEE", 0.0);
		vec_index++;
	   }
	}
	return 0;
}

int SBNsinglephoton::LoadSpectraApplyFullScaling(){

	//SBNspec temp_cv = *m_cv_spectrum;
	SBNspec temp_cv(tag+"_CV.SBNspec.root", this->xmlname, false);


	//m_scaled_spec_grid.clear();
	m_scaled_spec_grid.resize(m_total_gridpoints); 


	int ip_processd = 0; 
	if(m_bool_poly_grid){
	    std::cout<<"SBNsinglephoton::LoadSpectraApplyFullScaling\t|| Grab pre-scaled spectra and build the final spectra!" << std::endl; 
 	    //loop over polynomial grid
	    for(int i=0;i< m_poly_total_gridpoints; i++){
		std::vector<double> ipoint = m_vec_poly_grid[i];

		//first initilize SBNspec with pre-scaled root file
		std::ostringstream prescale_tag;
        	for(int ip=0; ip< ipoint.size(); ip++){
                	prescale_tag << "_" << std::fixed<< std::setprecision(3) << ipoint[ip];
                }	
		std::string full_filename = this->tag+"_PreScaled"+prescale_tag.str() + ".SBNspec.root";
		SBNspec temp_prescaled(full_filename, this->xmlname, false);	

		//loop over regular grid
		for(int j=0; j< m_num_total_gridpoints;j++){
			if(is_verbose && j%100==0 ) std::cout << "On Point " << i*m_num_total_gridpoints+j <<"/"<<m_total_gridpoints << std::endl;

			std::vector<double> jpoint = m_vec_grid[j];
		
			//m_scaled_spec_grid.push_back(new SBNspec((tag+"_CV.SBNspec.root").c_str(), this->xmlname, false));//start with genie CV
			m_scaled_spec_grid[ip_processd] = temp_cv; //start with genie CV
			for(int jp=0; jp< jpoint.size(); jp++){
				//m_scaled_spec_grid.back()->Scale(m_grid.f_dimensions[jp].f_name, jpoint[jp]); // scale corresponding subchannel
				m_scaled_spec_grid[ip_processd].Scale(m_grid.f_dimensions[jp].f_name, jpoint[jp]); // scale corresponding subchannel
			}
			
			//m_scaled_spec_grid.back()->Add(&temp_prescaled); //here we get the scaled spectra at final stage!!!
			m_scaled_spec_grid[ip_processd].Add(&temp_prescaled); //here we get the scaled spectra at final stage!!!

			// do not want to keep LEE in it
			if(tag == "NCpi0") m_scaled_spec_grid[ip_processd].Scale("NCDeltaRadOverlayLEE", 0.0);
			ip_processd++;

		}//end loop over regular grid

	     }// end loop over poly grid
	}
	else{
	    std::cout<<"SBNsinglephoton::LoadSpectraApplyFullScaling\t|| Scale the spectra for the whole grid " << std::endl;
	     for(int j=0; j< m_num_total_gridpoints;j++){
                        if(is_verbose && (j%1000 ==0)) std::cout << "On Point " << j <<"/"<<m_total_gridpoints << std::endl;

			std::vector<double> jpoint = m_vec_grid[j];

                        //m_scaled_spec_grid.push_back(new SBNspec((tag+"_CV.SBNspec.root").c_str(), this->xmlname, false));//start with genie CV
			m_scaled_spec_grid[ip_processd] = temp_cv; //start with genie CV
			for(int jp=0; jp< jpoint.size(); jp++){
                                //m_scaled_spec_grid.back()->Scale(m_grid.f_dimensions[jp].f_name, jpoint[jp]); // scale corresponding subchannel
                                m_scaled_spec_grid[ip_processd].Scale(m_grid.f_dimensions[jp].f_name, jpoint[jp]); // scale corresponding subchannel
                                //m_scaled_spec_grid[ip_processd].ScaleAll(jpoint[jp]); // scale corresponding subchannel
			}

	     		// do not want to keep LEE in it
			if(tag == "NCpi0") m_scaled_spec_grid[ip_processd].Scale("NCDeltaRadOverlayLEE", 0.0);
	
			ip_processd++;
	    }
	}

	return 0;
}



int SBNsinglephoton::CalcChiGridScanShapeOnlyFit(){
    if(!m_bool_cv_spectrum_loaded){
	std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| CV spec hasn't been loaded yet, load CV..." << std::endl;
	this->LoadCV();
    }	

    if(!m_bool_data_spectrum_loaded){
	std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| WARNING!! Data spec hasn't been loaded, will do a sensitivity study instead!" << std::endl;
	m_data_spectrum = new SBNspec();
	*m_data_spectrum = *m_cv_spectrum;
	m_data_spectrum->Scale("NCDeltaRadOverlayLEE", 0.0);
	m_data_spectrum->Scale("NCPi0Coh", 3.0);
	m_data_spectrum->Scale("NCPi0NotCoh", 0.8);
	SBNspec temp_data = *m_cv_spectrum;
	temp_data.Scale("NCDeltaRadOverlayLEE", 0.0);
	temp_data.CompareSBNspecs(m_data_spectrum, tag+"CVvsScaledCV");
	temp_data=*m_data_spectrum;
	this->PoissonFluctuation(m_data_spectrum);
	m_data_spectrum->WriteOut(tag+"_2g_fakedata");
	temp_data.CompareSBNspecs(m_data_spectrum, tag+"Before_AfterPoisson");

	temp_data = *m_cv_spectrum;
        temp_data.Scale("NCDeltaRadOverlayLEE", 0.0);
	temp_data.CompareSBNspecs(m_data_spectrum, tag+"CV_vs_FakeData");
    }else{
        m_cv_spectrum->CompareSBNspecs(m_data_spectrum, tag+"_CVvsData_NoErrorBar");
    }    


    if(m_scaled_spec_grid.size() != m_total_gridpoints){
	std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| ERROR!! # of scaled spectra: "<<m_scaled_spec_grid.size()<<" does NOT match grid size: " << m_total_gridpoints << " !!" <<std::endl;
	exit(EXIT_FAILURE);
    }

    //std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t" << __LINE__ << std::endl;

    SBNspec background_spectrum;
    SBNspec last_best_spectrum;
    double best_chi, last_best_chi;
    int best_point, last_best_point;
    std::vector<double> vec_chi, last_vec_chi;
    TFile* fout = new TFile(Form("%s_fit_output.root", tag.c_str()), "recreate");  //save matrix plot etc.

    m_full_but_genie_fractional_covariance_matrix->Write("frac_but_genie");
    m_genie_fractional_covariance_matrix->Write("frac_genie"); 

    m_chi = new SBNchi(this->xmlname);
    TMatrixT<double> full_systematic_covariance(num_bins_total, num_bins_total);
    TMatrixT<double> genie_systematic_matrix(num_bins_total, num_bins_total);
    TMatrixT<double> collapsed_full_systematic_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> inversed_total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);

    for(int n_iter =0; n_iter < m_max_number_iterations; n_iter ++){
	std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| On fit iteration "<< n_iter << std::endl;
	//reset everything at the biginning of each iteration
	best_chi = DBL_MAX;
	vec_chi.clear();


	if(n_iter == 0){
		last_best_spectrum = *m_cv_spectrum;
		last_best_spectrum.Scale("NCDeltaRadOverlayLEE", 0.0);
	}else{
		last_best_spectrum = m_scaled_spec_grid[last_best_point];		
	}

	background_spectrum = last_best_spectrum;


	//============================Calculate covariance matrix and its invert====================================/
	//full systematic covariance matrix, except genie uncertainty
	full_systematic_covariance = m_chi->FillSystMatrix(*m_full_but_genie_fractional_covariance_matrix, last_best_spectrum.full_vector);
	//calculate the shape only covariance matrix for genie uncertainty, to get rid of normalization uncertainty
	for(int i=0; i<m_grid.f_num_dimensions; i++){
	   //if(m_grid.f_dimensions[i].f_name == "NCDeltaRadOverlayLEE") continue;
	   SBNspec temp_comp = last_best_spectrum;
	   temp_comp.Keep( m_grid.f_dimensions[i].f_name, 1.0 );

	   genie_systematic_matrix = m_chi->FillSystMatrix(*m_genie_fractional_covariance_matrix, temp_comp.full_vector);
	   fout->cd();
	   genie_systematic_matrix.Write(Form("full_genie_%s_%d", m_grid.f_dimensions[i].f_name.c_str(), n_iter));	
	
	    //shape only matrix
	   //genie_systematic_matrix = m_chi->CalcShapeOnlyCovarianceMatrix(*m_genie_fractional_covariance_matrix, &temp_comp, &temp_comp);

	   //shape_mixed
	   genie_systematic_matrix = m_chi->CalcShapeMixedCovarianceMatrix(*m_genie_fractional_covariance_matrix, &temp_comp, &temp_comp);
	   genie_systematic_matrix.Write(Form("ShapeMixed_genie_%s_%d", m_grid.f_dimensions[i].f_name.c_str(), n_iter));	
	   full_systematic_covariance += genie_systematic_matrix;

	   background_spectrum.Scale(m_grid.f_dimensions[i].f_name, 0.0);   
	}
	//add genie uncertainties for other subchannels to total covariance matrix
	genie_systematic_matrix = m_chi->FillSystMatrix(*m_genie_fractional_covariance_matrix, background_spectrum.full_vector);
	full_systematic_covariance += genie_systematic_matrix;

	m_chi->CollapseModes(full_systematic_covariance, collapsed_full_systematic_matrix);

	//use CNP statistical covariance matrix
	total_covariance_matrix = m_chi->AddStatMatrixCNP(&collapsed_full_systematic_matrix, last_best_spectrum.collapsed_vector, m_data_spectrum->collapsed_vector);
	//std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t||check " << __LINE__ << std::endl;
	fout->cd();
	genie_systematic_matrix.Write(Form("full_genie_otherbkg_%d", n_iter));
	full_systematic_covariance.Write(Form("syst_uncollapsed_matrix_%d", n_iter));
	collapsed_full_systematic_matrix.Write(Form("syst_collapsed_matrix_%d", n_iter));
	total_covariance_matrix.Write(Form("total_collapsed_matrix_%d", n_iter));

	inversed_total_covariance_matrix= m_chi->InvertMatrix(total_covariance_matrix);


	if(is_verbose && m_bool_data_spectrum_loaded) last_best_spectrum.CompareSBNspecs(collapsed_full_systematic_matrix, m_data_spectrum, tag+"_Iter_"+std::to_string(n_iter));	
	//============================Done calculating covariance matrix ============================================/


	for(int i=0 ;i <m_total_gridpoints; i++){
	   if(i%1000 == 0) std::cout<< "On Point " << i << "/" << m_total_gridpoints << std::endl;
	   double temp_chi = m_chi->CalcChi(inversed_total_covariance_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector);
	   vec_chi.push_back(temp_chi);

	   if(temp_chi < best_chi){
		best_chi = temp_chi;
	        best_point = i;
	   }
	}

	if(is_verbose) std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| chi2 minimum :" << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;

	if(n_iter != 0){
	   if(fabs(best_chi - last_best_chi) < m_chi_min_convergance_tolerance){
		std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| chi2 has converged with best chi2 value " << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;
		break;
	   }
	}

	last_best_chi = best_chi;
	last_best_point = best_point;
	last_vec_chi = vec_chi;
    }//end loop for iteration

    fout->Close();
    if(is_verbose ){
       // collapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, m_scaled_spec_grid[best_point].full_vector, true);
        //SBNspec temp_best_spec = this->GeneratePointSpectra(best_point);
        SBNspec temp_best_spec = m_scaled_spec_grid[best_point];
	if(tag == "NCpi0") temp_best_spec.Scale("NCDeltaRadOverlayLEE", 0.0);
        temp_best_spec.CompareSBNspecs(collapsed_full_systematic_matrix, m_data_spectrum, tag+"_BFvsData");
    }
	
    //std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t||check " << __LINE__ << std::endl;
    m_map={{best_point, vec_chi}};
    this->PrintOutFitInfo(m_map, "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t||"+tag, true);
    this->WriteOutInfo(m_map);
    if(m_file_open) this->CloseFiles();
    return 0;
}



int SBNsinglephoton::CalcChiGridScan(){
    if(!m_bool_cv_spectrum_loaded){
	std::cout << "SBNsinglephoton::CalcChiGridScan\t|| CV spec hasn't been loaded yet, load CV..." << std::endl;
	this->LoadCV();
    }	

    if(!m_bool_data_spectrum_loaded){
	std::cout << "SBNsinglephoton::CalcChiGridScan\t|| WARNING!! Data spec hasn't been loaded, will do a sensitivity study instead!" << std::endl;
	m_data_spectrum = new SBNspec(tag+"_CV.SBNspec.root", xmlname, false);
	//*m_data_spectrum = *m_cv_spectrum;
	m_data_spectrum->Scale("NCPi0Coh", 3.0);
        m_data_spectrum->Scale("NCPi0NotCoh", 0.8);
	m_data_spectrum->Scale("NCDeltaRadOverlayLEE", 0.0);
	if(m_bool_modify_cv) m_data_spectrum->Scale("NCDeltaRadOverlayLEE", (m_cv_delta_scaling-1)*0.5);
	this->PoissonFluctuation(m_data_spectrum);
	m_data_spectrum->CollapseVector();
	m_data_spectrum->WriteOut("1g2g_fake_data");
	SBNspec temp_data = *m_cv_spectrum;
	temp_data.Scale("NCDeltaRadOverlayLEE", 0.0);
	temp_data.WriteOut(tag+"_CorrectedCV");
	temp_data.CompareSBNspecs(m_data_spectrum, tag+"_CorrectedCV_vs_Data");

    }else{
	m_cv_spectrum->CompareSBNspecs(m_data_spectrum, tag+"_CVvsData_NoErrorBar");
    }    


    if(m_scaled_spec_grid.size() != m_total_gridpoints){
	std::cout << "SBNsinglephoton::CalcChiGridScan\t|| ERROR!! # of scaled spectra: "<<m_scaled_spec_grid.size()<<" does NOT match grid size: " << m_total_gridpoints << " !!" <<std::endl;
	exit(EXIT_FAILURE);
    }


    SBNspec last_best_spectrum;
    double best_chi, last_best_chi;
    int best_point, last_best_point;
    std::vector<double> vec_chi, last_vec_chi;
    TFile* fout = new TFile(Form("%s_fit_output.root", tag.c_str()), "recreate");  //save matrix plot etc.

    m_full_fractional_covariance_matrix->Write("full_fractional_matrix");
    //std::cout << "SBNsinglephoton::CalcChiGridScan\t||check " << __LINE__ << std::endl;

    m_chi = new SBNchi(this->xmlname);
    TMatrixT<double> collapsed_full_systematic_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> inversed_total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);

    TMatrixT<double>* collapse_frac;

    for(int n_iter =0; n_iter < m_max_number_iterations; n_iter ++){
	std::cout << "SBNsinglephoton::CalcChiGridScan\t|| On fit iteration "<< n_iter << std::endl;
	//reset everything at the biginning of each iteration
	best_chi = DBL_MAX;
	vec_chi.clear();


	if(n_iter == 0){
		last_best_spectrum = *m_cv_spectrum;
		last_best_spectrum.Scale("NCDeltaRadOverlayLEE", 0.0);
	}else{
		last_best_spectrum = m_scaled_spec_grid[last_best_point];		
	}


	//============================Calculate covariance matrix and its invert====================================/
	//full systematic covariance matrix, except genie uncertainty

	//for(int i=0; i<last_best_spectrum.full_vector.size(); i++)
	//	std::cout << last_best_spectrum.full_vector[i] << ",";
	//std::cout << std::endl;

	//SBNspec temp_signal = last_best_spectrum; temp_signal.Keep("NCDeltaRadOverlaySM", 1.0);
	//SBNspec temp_other = last_best_spectrum; temp_other.Scale("NCDeltaRadOverlaySM", 0.0);
	//collapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, temp_signal.full_vector, true) + m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, temp_other.full_vector, true);

	collapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, last_best_spectrum.full_vector, true);
	collapse_frac = (TMatrixT<double>*)collapsed_full_systematic_matrix.Clone();
	for(int i=0; i<num_bins_total_compressed ;i ++){
	   for(int j=0 ; j< num_bins_total_compressed; j++)
		(*collapse_frac)(i,j) /= last_best_spectrum.collapsed_vector[i]* last_best_spectrum.collapsed_vector[j];
	}
	//use CNP statistical covariance matrix
	//total_covariance_matrix = m_chi->AddStatMatrix(&collapsed_full_systematic_matrix, m_data_spectrum->collapsed_vector);
	total_covariance_matrix = m_chi->AddStatMatrixCNP(&collapsed_full_systematic_matrix, last_best_spectrum.collapsed_vector, m_data_spectrum->collapsed_vector);

	inversed_total_covariance_matrix= m_chi->InvertMatrix(total_covariance_matrix);

	fout->cd();
	collapsed_full_systematic_matrix.Write(Form("syst_collapsed_matrix_%d", n_iter));
	total_covariance_matrix.Write(Form("total_collapsed_matrix_%d", n_iter));
	collapse_frac->Write(Form("collapsed_fractional_matrix_%d", n_iter));
	inversed_total_covariance_matrix.Write(Form("inversed_matrix_%d", n_iter));

	if(is_verbose && m_bool_data_spectrum_loaded) last_best_spectrum.CompareSBNspecs(collapsed_full_systematic_matrix, m_data_spectrum, tag+"_Iter_"+std::to_string(n_iter));	
	//============================Done calculating covariance matrix ============================================/


	for(int i=0 ;i <m_total_gridpoints; i++){
	   if(i%1000 ==0) std::cout<< "On Point " << i << "/" << m_total_gridpoints << std::endl;
	   double temp_chi;
	   //if(i == 400 || i== 950){ 
	   //	std::cout << " On point " << i << std::endl;
	   //	temp_chi = m_chi->CalcChi(inversed_total_covariance_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector, true);
	   //}
	   //else temp_chi = m_chi->CalcChi(inversed_total_covariance_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector, false);
	   temp_chi = m_chi->CalcChi(inversed_total_covariance_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector, false);
	   vec_chi.push_back(temp_chi);

	   if(temp_chi < best_chi){
		best_chi = temp_chi;
	        best_point = i;
	   }
	}

	if(is_verbose) std::cout << "SBNsinglephoton::CalcChiGridScan\t|| chi2 minimum :" << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;

	if(n_iter != 0){
	   if(fabs(best_chi - last_best_chi) < m_chi_min_convergance_tolerance){
		std::cout << "SBNsinglephoton::CalcChiGridScan\t|| chi2 has converged with best chi2 value " << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;
		break;
	   }
	}

	last_best_chi = best_chi;
	last_best_point = best_point;
	last_vec_chi = vec_chi;
    }//end loop for iteration

    fout->Close();

    //best-fit vs data comparison	
    //if(is_verbose && m_bool_data_spectrum_loaded){
    if(is_verbose ){
        collapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, m_scaled_spec_grid[best_point].full_vector, true);
	SBNspec temp_best_spec = this->GeneratePointSpectra(best_point);
	if(tag == "NCpi0") temp_best_spec.Scale("NCDeltaRadOverlayLEE", 0.0);
        temp_best_spec.CompareSBNspecs(collapsed_full_systematic_matrix, m_data_spectrum, tag+"_BFvsData");
    }

   std::cout << " BEST-FIT \t\t\t Data \t\t\t CV " << std::endl;
   for(int i=0 ; i< m_data_spectrum->collapsed_vector.size(); i++){
	std::cout <<i << ": " << m_scaled_spec_grid[best_point].collapsed_vector[i] << "\t\t\t " << m_data_spectrum->collapsed_vector[i] << " \t\t\t " << m_cv_spectrum->collapsed_vector[i] << std::endl;
   }
    m_map={{best_point, vec_chi}};
    this->PrintOutFitInfo(m_map, "SBNsinglephoton::CalcChiGridScan\t|| "+tag, true);
    this->WriteOutInfo(m_map);
    if(m_file_open) this->CloseFiles();
    return 0;
}


int SBNsinglephoton::LoadCV(){
    if(is_verbose) std::cout << "SBNsinglephoton::LoadCV\t|| Setup CV spectrum" << std::endl;
    m_cv_spectrum = new SBNspec(tag+"_CV.SBNspec.root", xmlname, false);
    m_cv_spectrum->CollapseVector();
//    m_cv_spectrum->Scale("NCPi0Coh", 1.25);
//    m_cv_spectrum->Scale("NCPi0NotCoh", 0.9);
//    m_cv_spectrum->CollapseVector();
    m_bool_cv_spectrum_loaded = true;


    return 0;
}

int SBNsinglephoton::LoadData(std::string filename){
    if(is_verbose) std::cout << "SBNsinglephoton::LoadData\t|| Load data spectrum from file: " << filename << std::endl;
    m_data_spectrum = new SBNspec(filename, xmlname, false);
    m_bool_data_spectrum_loaded = true;
    m_data_spectrum->CollapseVector(); 
    return 0;
}


int SBNsinglephoton::SetFullFractionalCovarianceMatrix(std::string filename, std::string matrix_name){

	if(is_verbose) std::cout << "SBNsinglephoton::SetFullFractionalCovarianceMatrix||\tOpen total frac covariance matrix file: " << filename << std::endl;
	TFile* f_syst = new TFile(filename.c_str(), "read");
	m_full_fractional_covariance_matrix = (TMatrixT<double>*)f_syst->Get(matrix_name.c_str());

	if(m_full_fractional_covariance_matrix->GetNcols() != m_full_fractional_covariance_matrix->GetNrows()){
		std::cout << "SBNsinglephoton::SetFullFractionalCovarianceMatrix\t|| Matrix provided is not sysmetric" << std::endl;
		exit(EXIT_FAILURE);		
	}

	this->RemoveNan(m_full_fractional_covariance_matrix);


/*	//get a subset of matrix
	TMatrixT<double> temp_matrix = *m_full_fractional_covariance_matrix;
	m_full_fractional_covariance_matrix->ResizeTo(num_bins_total, num_bins_total);
	// *m_full_fractional_covariance_matrix  = temp_matrix.GetSub(90, 329, 90, 329, "S");
	// *m_full_fractional_covariance_matrix  = temp_matrix.GetSub(0, 89, 0, 89, "S");
	*m_full_fractional_covariance_matrix  = temp_matrix.GetSub(0, 29, 0, 29, "S");
	
	//get submatrices from a big covariance matrix
	//for 1g1p+2g1p
	m_full_fractional_covariance_matrix->SetSub(0, 0, temp_matrix.GetSub(0, 29, 0, 29, "S"));
	m_full_fractional_covariance_matrix->SetSub(0, 30, temp_matrix.GetSub(0, 29, 90, 209, "S"));
	m_full_fractional_covariance_matrix->SetSub(30, 0, temp_matrix.GetSub(90, 209, 0, 29, "S"));
	m_full_fractional_covariance_matrix->SetSub(30, 30, temp_matrix.GetSub(90, 209, 90, 209, "S"));
	//for 1g0p+2g0p
	m_full_fractional_covariance_matrix->SetSub(0, 0, temp_matrix.GetSub(30, 89, 30, 89, "S"));
	m_full_fractional_covariance_matrix->SetSub(0, 60, temp_matrix.GetSub(30, 89, 210, 329, "S"));
	m_full_fractional_covariance_matrix->SetSub(60, 0, temp_matrix.GetSub(210, 329, 30, 89, "S"));
	m_full_fractional_covariance_matrix->SetSub(60, 60, temp_matrix.GetSub(210, 329, 210, 329, "S"));
*/	std::cout << "total bins " << num_bins_total << ", matrix size " << m_full_fractional_covariance_matrix->GetNrows() << std::endl;
	this->RemoveNan(m_full_fractional_covariance_matrix);

	if(is_verbose)std::cout << "SBNsinglephoton::SetFullFractionalCovarianceMatrix\t|| matrix size: " << m_full_fractional_covariance_matrix->GetNcols()<< std::endl;;

/*
	TMatrixT<double>* temp_genie_XS = (TMatrixT<double>*)f_syst->Get("individualDir/All_UBGenie_frac_covariance");
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/AxFFCCQEshape_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/DecayAngMEC_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/NormCCCOH_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/NormNCCOH_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/RPA_CCQE_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/Theta_Delta2Npi_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/VecFFCCQEshape_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/XSecShape_CCMEC_UBGenie_frac_covariance"));
        this->RemoveNan(temp_genie_XS);

	*m_full_fractional_covariance_matrix = *temp_genie_XS;

	for(int i=0; i<m_full_fractional_covariance_matrix->GetNcols(); i++){
	    for(int j=0; j<m_full_fractional_covariance_matrix->GetNcols(); j++){
		//if(std::isnan((*m_full_fractional_covariance_matrix)(i,j))) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
		if( i>89 && j <= 89) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
		if( i<=89 && j > 89) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
		//if( (i== 2 || i==8) && j == 329) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
		//if( (j== 2 || j==8) && i == 329) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
	    }	
	}
*/
	f_syst->Close();
	return 0;
}

int SBNsinglephoton::SetGenieFractionalCovarianceMatrix(std::string filename){

	//check with Gray, are matrices added directly before checking the nan's?
	if(is_verbose) std::cout <<"SBNsinglephoton::SetGenieFractionalCovarianceMatrix||\tOpen genie matrix file: " << filename << std::endl;
	TFile* f_syst = new TFile(filename.c_str(), "read");
	//TFile* f_syst = new TFile(filename.c_str(), "UPDATE");
	//m_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
	m_genie_fractional_covariance_matrix  = (TMatrixT<double>*)f_syst->Get("individualDir/All_UBGenie_frac_covariance");
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/AxFFCCQEshape_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/DecayAngMEC_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/NormCCCOH_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/NormNCCOH_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/RPA_CCQE_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/Theta_Delta2Npi_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/VecFFCCQEshape_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/XSecShape_CCMEC_UBGenie_frac_covariance"));
	this->RemoveNan(m_genie_fractional_covariance_matrix);


/*
	SBNspec temp_cv = *m_cv_spectrum;
	temp_cv.Keep("NCPi0Coh", 1.0);
	temp_cv.CalcFullVector();
	std::vector<double> temp_full = temp_cv.full_vector;


	TMatrixT<double> M_cohOnly = *m_genie_fractional_covariance_matrix;
	for(int i=0; i<num_bins_total; i++){
	   for(int j=0; j< num_bins_total; j++){
	     if(temp_full[i]==0 && temp_full[j]==0 ) M_cohOnly(i,j) = 0;
	   }
	}

	TMatrixT<double>* M_updatedCohOnly = (TMatrixT<double>*)f_syst->Get("individualDir/updated_allUBGenie_frac_covariance");
	for(int i=0; i<num_bins_total; i++){
           for(int j=0; j< num_bins_total; j++){
             if(temp_full[i]==0 && temp_full[j]==0 ) (*M_updatedCohOnly)(i,j) = 0;
           }
        }


	TCanvas c("c", "c",100, 100);
	gStyle->SetOptStat(0);
	M_cohOnly.Draw("colz");
	c.Update();
	c.SaveAs("Original_Genie_Matrix_COHonly.pdf", "pdf");
	c.Clear();
	c.cd();
	M_updatedCohOnly->Draw("colz");
	c.Update();
	c.SaveAs("Updated_Genie_Matrix_COHonly.pdf", "pdf");
*/

/*	//update the genie uncertainty of COH pi0's to be 1.
	TMatrixT<double> temp_matrix = *m_genie_fractional_covariance_matrix;
	SBNspec temp_cv = *m_cv_spectrum;
	temp_cv.Keep("NCPi0Coh", 1.0);
	temp_cv.CalcFullVector();
	std::vector<double> temp_full = temp_cv.full_vector;

	for(int i=0; i<temp_full.size(); i++){
	   for(int j=0;j<temp_full.size(); j++){
		if(temp_full[i]!=0 && temp_full[j]!=0) (*m_genie_fractional_covariance_matrix)(i,j)=1.0;
	   }
	}
	std::cout << "check " << __LINE__ << std::endl;

	TMatrixT<double>* m_fluxxs = (TMatrixT<double>*)f_syst->Get("frac_covariance");
	this->RemoveNan(m_fluxxs);
	*m_fluxxs +=  *m_genie_fractional_covariance_matrix - temp_matrix;
	f_syst->cd();
	m_fluxxs->Write("updated_frac_covariance", TObject::kWriteDelete);
	TDirectory *individualDir = f_syst->GetDirectory("individualDir");
	individualDir->cd();
	m_genie_fractional_covariance_matrix->Write("updated_allUBGenie_frac_covariance", TObject::kWriteDelete);
	std::cout << "check " << __LINE__ << std::endl;

	for(int i=0; i<num_bins_total; i++){
	    for(int j=0; j< num_bins_total; j++){
		if(temp_matrix(i,j) != (*m_genie_fractional_covariance_matrix)(i,j)){
		    if(temp_full[i]==0 || temp_full[j] ==0)
			std::cout << "We got trouble" << std::endl;
		}
	    }
	}


	 TFile* f_all = new TFile("/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/FakeData/Poisson_10/updated_1g1p_1g0p_2g1p_2g0p_Combined_FluxXSDet_v11_Momentum_fracfixed.SBNcovar.root", "UPDATE");
         TMatrixT<double>* mall = (TMatrixT<double>*)f_all->Get("frac_covariance");
	this->RemoveNan(mall);
	TMatrixT<double> mall_temp  = *mall;
	mall_temp += *m_genie_fractional_covariance_matrix - temp_matrix;
	f_all->cd();
	mall_temp.Write("updated_frac_covariance",TObject::kWriteDelete);
	f_all->Close();
	std::cout << "check " << __LINE__ << std::endl;
*/


/*	//set genie matrix to be the same as total covariance matrix
	m_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
	*m_genie_fractional_covariance_matrix = *m_full_fractional_covariance_matrix;
*/


	//get a sub of matrix
/*	TMatrixT<double> temp_matrix = *m_genie_fractional_covariance_matrix;
	m_genie_fractional_covariance_matrix->ResizeTo(num_bins_total, num_bins_total);
	*m_genie_fractional_covariance_matrix  = temp_matrix.GetSub(90, 329, 90, 329, "S");
*/
	std::cout << "total bins " << num_bins_total << ", matrix size " << m_genie_fractional_covariance_matrix->GetNrows() << std::endl;

//	this->RemoveNan(m_genie_fractional_covariance_matrix);

	TMatrixT<double>* temp = (TMatrixT<double>*)f_syst->Get("frac_covariance");
	this->RemoveNan(temp);
	*m_full_fractional_covariance_matrix += *temp;	

	delete temp;

	f_syst->Close();
	
	return 0;
}


int SBNsinglephoton::CalcFullButGenieFractionalCovarMatrix(){

	if(m_full_fractional_covariance_matrix==NULL || m_genie_fractional_covariance_matrix==NULL){
	   std::cout<< "SBNsinglephoton::CalcFullButGenieFractionalCovarMatrix\t|| Either full fractional covar matrix or genie fractional covariance matrix has NOT been setup yet."<< std::endl;
	   exit(EXIT_FAILURE); 
	}

	if(is_verbose) std::cout << "SBNsinglephoton::CalcFullButGenieFractionalCovarMatrix\t|| as the name says.." << std::endl;
	m_full_but_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
	*m_full_but_genie_fractional_covariance_matrix  = (*m_full_fractional_covariance_matrix) - (*m_genie_fractional_covariance_matrix);

	return 0;
}


double SBNsinglephoton::ScaleFactor(double E, double factor, std::vector<double>& vec){
	double scale = factor;
	switch(vec.size()){
	   case 0:
		//std::cout<< "SBNsinglephoton::ScaleFactor\t|| No energy/momentum dependent scaling applied" << std::endl;
		break;
	   case 1:
		//std::cout << "SBNsinglephoton::ScaleFactor\t|| Applying Linear energy/momentum dependent scaling!";
		scale += vec[0]*E;
		//std::cout << " " << scale << std::endl;
		break;
	   case 2:
		//std::cout << "SBNsinglephoton::ScaleFactor\t|| Applying 2nd order polynomial energy/momentum dependent scaling!" << std::endl;
		scale += vec[0]*E+vec[1]*E*E;
		break;
	   default: scale += 0;
	}

	//if(scale < 0) scale = 1.0;  //don't allow negative weight
	if(scale < 0) scale = 0;  //don't allow negative weight
	return scale;

}


int SBNsinglephoton::RemoveNan(TMatrixT<double>* M){

	int row_size = M->GetNrows();
	int col_size = M->GetNcols();

	if(is_verbose) std::cout<< "SBNsinglephoton::RemoveNan\t|| Remove the nan's from matrix " << M->GetName() << std::endl;
        for(int i=0; i< row_size; i++){
            for(int j=0; j< col_size; j++){
                if(std::isnan((*M)(i,j))) (*M)(i,j) =0.0;
            }
        }

	return 0;	
}


int SBNsinglephoton::GrabFitMap(){
	if(is_verbose) std::cout << "SBNsinglephoton::GrabFitMap\t||\tOpen file: "<< tag<<"_fit_output.root, and grab chi and BF index info"<< std::endl;
	fin = new TFile(Form("%s_fit_output.root", tag.c_str()), "UPDATE");
	TVectorD* v = (TVectorD*)fin->Get("bf_index");	
	int bf_index = (int)(*v)[0];
	std::vector<double>* p_vec_chi;
	fin->GetObject("vector_chi", p_vec_chi);

	m_map={{bf_index, *p_vec_chi}};		

	delete v;
	delete p_vec_chi;
	return 0;
}


int SBNsinglephoton::SaveHistogram(){	
        return this->SaveHistogram(m_map);
}

int SBNsinglephoton::SaveHistogram(std::map<int, std::vector<double>>& inmap){


    if(inmap.empty()){
	std::cout << "SBNsinglephoton::SaveHistogram\t|| map is empty!!" << std::endl;
	exit(EXIT_FAILURE);	
    }else{
        std::map<int, std::vector<double>> map_chi = inmap;
        std::map<int, std::vector<double>>::iterator itmap = map_chi.begin();
	int best_point = itmap->first;
	std::vector<double> vec_chi = itmap->second;
	double chi_min = *std::min_element(vec_chi.begin(), vec_chi.end());
	std::for_each(vec_chi.begin(), vec_chi.end(), [&chi_min](double& d){d -= chi_min;});  //get the delta_chi vector


	std::map<std::string, std::string> title_map={{"NCDeltaRadOverlaySM", "Factor for NC #Delta Radiative x_{#Delta}"},
						      {"NCDeltaRadOverlayLEE", "Factor for NC #Delta Radiative x_{#Delta}"},
						      {"NCPi0Coh", "Factor for NC 1 #pi^{0} Coherent"},
						      {"NCPi0NotCoh", "Factor for NC 1 #pi^{0} Non-Coherent"},
						      {"All", "Flat normalization factor for all"}};

        int m_poly_best_index = best_point/m_num_total_gridpoints;
        int m_best_index = best_point%m_num_total_gridpoints;

	  
	if(m_grid.f_num_dimensions == 2){
	
	   if(!m_bool_poly_grid){
	    	if(is_verbose) std::cout<< "SBNsinglephoton::SaveHistogram\t|| Case: NCpi0 normalization fit, no energy/momentum dependent scaling!" << std::endl;
		auto grid_x = m_grid.f_dimensions.at(0);
		auto grid_y = m_grid.f_dimensions.at(1);

		TH2D h_chi_surface = this->Do2DInterpolation(m_interpolation_number, grid_y.f_points, grid_x.f_points,vec_chi, tag);
		h_chi_surface.SetName("h_chi_interpolated_surface");
		h_chi_surface.SetTitle(Form("#Delta#chi^{2} surface; %s;%s", title_map[grid_y.f_name].c_str(), title_map[grid_x.f_name].c_str()));
		h_chi_surface.Write();
	        TH2D* hr=(TH2D*)h_chi_surface.Clone(); 

		//draw contours
		std::vector<TGraph*> contour_graph = this->FindContour(h_chi_surface, 3, tag);  //want 3 contour
	        //TCanvas* c_canvas=(TCanvas*)this->DrawContour(hr, contour_graph, m_vec_grid[m_best_index]).Clone();	
	        this->DrawContour(hr, contour_graph, m_vec_grid[m_best_index]);	
		
	   }else{
		if(is_verbose) std::cout<< "SBNsinglephoton::SaveHistogram\t|| Case: NCpi0 normalization fit, with energy/momentum dependent scaling to " << m_poly_grid.f_num_dimensions << "nd order!" << std::endl;
		 auto grid_1order = m_poly_grid.f_dimensions.at(0);

		 //grab the flat grid associated with poly_grid
		 NGridDimension grid_flat = grid_1order;
		 NGridDimension grid_fother = grid_1order;
		 int period = 0;  // very important variable!! 
		 for(auto &grid:m_grid.f_dimensions){
			period *= grid.f_N;
			if(grid.f_name == grid_1order.f_name){
			    grid_flat = grid;
			    period = 1;
			}else{
			    grid_fother = grid;
			}
		 }		
		 TH2D h_mchi_poly("h_mchi_poly", Form("marginalized #Delta#chi^{2}; %s flat; %s 1st order", title_map[grid_flat.f_name].c_str(), title_map[grid_flat.f_name].c_str()), grid_flat.f_N, grid_flat.f_min, grid_flat.f_max, grid_1order.f_N, grid_1order.f_min, grid_1order.f_max);
		 TH2D h_mchi_flat("h_mchi_flat", Form("marginalized #Delta#chi^{2}; %s; %s", title_map[grid_fother.f_name].c_str(), title_map[grid_flat.f_name].c_str()), grid_fother.f_N, grid_fother.f_min, grid_fother.f_max, grid_flat.f_N, grid_flat.f_min, grid_flat.f_max);

		 TH2D h_mchi_mix("h_mchi_mix", Form("marginalized #Delta#chi^{2}; %s; %s 1st order", title_map[grid_fother.f_name].c_str(), title_map[grid_1order.f_name].c_str()), grid_fother.f_N, grid_fother.f_min, grid_fother.f_max, grid_1order.f_N, grid_1order.f_min, grid_1order.f_max);
		 std::vector<TH1D*> vh_mchi(3);
		 vh_mchi[0] = new TH1D(Form("h_mchi_flat_%s", (grid_flat.f_name).c_str()), Form("marginalized #Delta#chi^{2}; %s;marginalized #Delta#chi^{2}",title_map[grid_flat.f_name].c_str()), grid_flat.f_N, grid_flat.f_min, grid_flat.f_max);
		 vh_mchi[1] = new TH1D(Form("h_mchi_poly_%s", (grid_1order.f_name).c_str()), Form("marginalized #Delta#chi^{2}; %s 1st order;marginalized #Delta#chi^{2}",title_map[grid_1order.f_name].c_str()), grid_1order.f_N, grid_1order.f_min, grid_1order.f_max);
		 vh_mchi[2] = new TH1D(Form("h_mchi_%s", (grid_fother.f_name).c_str()), Form("marginalized #Delta#chi^{2}; %s; marginalized #Delta#chi^{2}",title_map[grid_fother.f_name].c_str()),grid_fother.f_N, grid_fother.f_min, grid_fother.f_max);

		 for(int i=0;i<grid_flat.f_N; i++){
			vh_mchi[0]->SetBinContent(i+1, DBL_MAX);
			for(int j=0;j<grid_1order.f_N; j++)
			    h_mchi_poly.SetBinContent(i+1, j+1, DBL_MAX);
			for(int j=0;j<grid_fother.f_N; j++)
			    h_mchi_flat.SetBinContent(j+1, i+1, DBL_MAX);
		 } 	
		 for(int i=0;i<grid_fother.f_N; i++){
			vh_mchi[2]->SetBinContent(i+1,DBL_MAX);
			for(int j=0;j<grid_1order.f_N; j++)
			    h_mchi_mix.SetBinContent(i+1, j+1, DBL_MAX);
		 } 
		for(int i=0; i<grid_1order.f_N; i++) vh_mchi[1]->SetBinContent(i+1, DBL_MAX);
	
		 for(int i=0; i<vec_chi.size(); i++){
		     int temp_poly_index = i/m_num_total_gridpoints;
        	     int temp_mgrid_index = i%m_num_total_gridpoints; 
		     int flat_index = (temp_mgrid_index%(period* grid_flat.f_N))/period; //get its index in flat dimension for current point
			// needs more work
		     int fother_index = temp_mgrid_index%grid_fother.f_N; //index for another flat dimension
		     if(period == 1) fother_index = temp_mgrid_index/grid_flat.f_N;
		     if(vec_chi[i] < h_mchi_poly.GetBinContent(flat_index+1, temp_poly_index+1)) h_mchi_poly.SetBinContent(flat_index+1, temp_poly_index+1, vec_chi[i]);
		     if(vec_chi[i] < h_mchi_flat.GetBinContent(fother_index+1, flat_index+1)) h_mchi_flat.SetBinContent(fother_index+1, flat_index+1, vec_chi[i]);
		     if(vec_chi[i] < h_mchi_mix.GetBinContent(fother_index+1, temp_poly_index+1)) h_mchi_mix.SetBinContent(fother_index+1, temp_poly_index+1, vec_chi[i]);
		     if(vec_chi[i] < vh_mchi[0]->GetBinContent(flat_index+1)) vh_mchi[0]->SetBinContent(flat_index+1, vec_chi[i]);
		     if(vec_chi[i] < vh_mchi[1]->GetBinContent(temp_poly_index+1)) vh_mchi[1]->SetBinContent(temp_poly_index+1, vec_chi[i]);
		     if(vec_chi[i] < vh_mchi[2]->GetBinContent(fother_index+1)) vh_mchi[2]->SetBinContent(fother_index+1, vec_chi[i]);
		 }	 
		fin->cd(); 
		h_mchi_poly.Write();
		h_mchi_flat.Write();
		h_mchi_mix.Write();

		for(auto &h:vh_mchi){
			h->Write();
			TCanvas* c=new TCanvas("c", "c");
			gStyle->SetOptStat(0);
			h->Draw("hist");
			TLine line(h->GetXaxis()->GetXmin(), 1.0, h->GetXaxis()->GetXmax(), 1.0);
			TLine line90(h->GetXaxis()->GetXmin(), 2.71, h->GetXaxis()->GetXmax(), 2.71);
			line.SetLineColor(8);
			line90.SetLineColor(42);
			line.Draw("same");
			line90.Draw("same");
			c->Update();
			c->SaveAs((tag+"_"+h->GetName()+".pdf").c_str(), "pdf");

		}
		//save marginalized chi
		std::vector<double> marginalized_chi_poly, marginalized_chi_mix;
		for(int i=0; i< grid_1order.f_N;i++){
		    for(int j=0; j< grid_flat.f_N; j++)
			marginalized_chi_poly.push_back(h_mchi_poly.GetBinContent(j+1, i+1));
		    for(int j=0; j< grid_fother.f_N; j++)
			marginalized_chi_mix.push_back(h_mchi_mix.GetBinContent(j+1, i+1));
		}


		//save marginalized chi
		std::vector<double> marginalized_chi_flat;
		for(int i=0; i< grid_flat.f_N;i++)
		    for(int j=0; j< grid_fother.f_N; j++)
			marginalized_chi_flat.push_back(h_mchi_flat.GetBinContent(j+1, i+1));

		//draw the contour
		TH2D h_mchi_flatinter = this->Do2DInterpolation(m_interpolation_number, grid_fother.f_points, grid_flat.f_points, marginalized_chi_flat, tag+"_Flat");
		std::vector<TGraph*> vg_mchi_flat_contour = this->FindContour(h_mchi_flatinter, 3, tag+"_Flat");
		h_mchi_flat.GetXaxis()->SetRangeUser(h_mchi_flatinter.GetXaxis()->GetXmin(), h_mchi_flatinter.GetXaxis()->GetXmax());
                h_mchi_flat.GetYaxis()->SetRangeUser(h_mchi_flatinter.GetYaxis()->GetXmin(), h_mchi_flatinter.GetYaxis()->GetXmax());
		this->DrawContour(&h_mchi_flat, vg_mchi_flat_contour, tag+"_Flat",std::vector<double>{});

		//draw the contour
		TH2D h_mchi_polyinter = this->Do2DInterpolation(m_interpolation_number, grid_flat.f_points, grid_1order.f_points, marginalized_chi_poly, tag+"_Poly");
		std::vector<TGraph*> vg_mchi_contour = this->FindContour(h_mchi_polyinter, 3, tag+"_Poly");
		h_mchi_poly.GetXaxis()->SetRangeUser(h_mchi_polyinter.GetXaxis()->GetXmin(), h_mchi_polyinter.GetXaxis()->GetXmax());
                h_mchi_poly.GetYaxis()->SetRangeUser(h_mchi_polyinter.GetYaxis()->GetXmin(), h_mchi_polyinter.GetYaxis()->GetXmax());
		this->DrawContour(&h_mchi_poly, vg_mchi_contour,tag+"_Poly", std::vector<double>{});

		TH2D h_mchi_mixinter = this->Do2DInterpolation(m_interpolation_number, grid_fother.f_points, grid_1order.f_points, marginalized_chi_mix, tag+"_Mix");
		std::vector<TGraph*> vg_mchi_mix_contour = this->FindContour(h_mchi_mixinter, 3, tag+"_Mix");
		h_mchi_mix.GetXaxis()->SetRangeUser(h_mchi_mixinter.GetXaxis()->GetXmin(), h_mchi_mixinter.GetXaxis()->GetXmax());
                h_mchi_mix.GetYaxis()->SetRangeUser(h_mchi_mixinter.GetYaxis()->GetXmin(), h_mchi_mixinter.GetYaxis()->GetXmax());
		this->DrawContour(&h_mchi_mix, vg_mchi_mix_contour,tag+"_Mix", std::vector<double>{});


	   } //end of if_poly_grid loop
	} //end of 2 dimension case
	else if(m_grid.f_num_dimensions == 1){
	   TH1D* h_dchi=nullptr;
	   NGridDimension xgrid = m_grid.f_dimensions.at(0);
	   if(xgrid.f_name == "NCDeltaRadOverlaySM") h_dchi = new TH1D("h_delta_chi", Form("#Delta#chi^{2} distribution;%s; #Delta#chi^{2} ",title_map[xgrid.f_name].c_str()), m_grid.f_num_total_points, xgrid.f_min, xgrid.f_max);
	   else if(xgrid.f_name == "NCDeltaRadOverlayLEE" ) h_dchi = new TH1D("h_delta_chi", Form("#Delta#chi^{2} distribution;%s; #Delta#chi^{2} ",title_map[xgrid.f_name].c_str()), m_grid.f_num_total_points, (xgrid.f_min)*2+1, (xgrid.f_max)*2+1);
	   for(int i=0 ;i< vec_chi.size(); i++){
		std::vector<double> ipoint = m_vec_grid[i];
		//h_dchi->Fill(ipoint[0], vec_chi[i]);
		h_dchi->SetBinContent(i+1, vec_chi[i]);
	   }

	   h_dchi->Write();
	   TCanvas c("c_chi_delta", "c_chi_delta");
                h_dchi->Draw("hist");
                TLine line(h_dchi->GetXaxis()->GetXmin(), 1.0, h_dchi->GetXaxis()->GetXmax(), 1.0);
                TLine line90(h_dchi->GetXaxis()->GetXmin(), 2.71, h_dchi->GetXaxis()->GetXmax(), 2.71);
                line.SetLineColor(8);
                line90.SetLineColor(42);
                line.Draw("same");
                line90.Draw("same");
                c.Update();
		c.SaveAs((tag+"_pretty_chi.pdf").c_str(), "pdf");
                c.Write();
	}//end of 1 dimension case
	else{

	   if(!m_bool_poly_grid){
		if(is_verbose) std::cout<< "SBNsinglephoton::SaveHistogram\t|| Case: NC delta and pi0 combined fit, no energy/momentum dependent scaling!" << std::endl;
		auto grid_x = m_grid.f_dimensions.at(0);
		auto grid_y = m_grid.f_dimensions.at(1);
		auto grid_z = m_grid.f_dimensions.at(2);   //assume grid_Z is the grid for NCDeltaRadOverlayLEE here.
		std::vector<double> temp_best_point = m_vec_grid[m_best_index];

		//marginalize over 1 parameter	
		TH2D* h_mchi2_xy = new TH2D("h_mchi2_xy", Form("marginalized #Delta#chi^{2} surface; %s;%s", title_map[grid_x.f_name].c_str(), title_map[grid_y.f_name].c_str()), grid_x.f_N, grid_x.f_min, grid_x.f_max, grid_y.f_N, grid_y.f_min, grid_y.f_max);
		TH2D* h_mchi2_yz = new TH2D("h_mchi2_yz", Form("marginalized #Delta#chi^{2} surface; %s;%s", title_map[grid_y.f_name].c_str(), title_map[grid_z.f_name].c_str()), grid_y.f_N, grid_y.f_min, grid_y.f_max, grid_z.f_N, (grid_z.f_min)*2+1, (grid_z.f_max)*2+1);
		TH2D* h_mchi2_xz = new TH2D("h_mchi2_xz", Form("marginalized #Delta#chi^{2} surface; %s;%s", title_map[grid_x.f_name].c_str(), title_map[grid_z.f_name].c_str()), grid_x.f_N, grid_x.f_min, grid_x.f_max, grid_z.f_N, (grid_z.f_min)*2+1, (grid_z.f_max)*2+1);
		//global minimum
		TH2D* h_gchi2_xy = new TH2D("h_gchi2_xy", Form("h_gchi2_xy; %s;%s", title_map[grid_x.f_name].c_str(), title_map[grid_y.f_name].c_str()), grid_x.f_N, grid_x.f_min, grid_x.f_max, grid_y.f_N, grid_y.f_min, grid_y.f_max);
		TH2D* h_gchi2_yz = new TH2D("h_gchi2_yz", Form("h_gchi2_yz; %s;%s", title_map[grid_y.f_name].c_str(), title_map[grid_z.f_name].c_str()), grid_y.f_N, grid_y.f_min, grid_y.f_max, grid_z.f_N, (grid_z.f_min)*2+1, (grid_z.f_max)*2+1);
		TH2D* h_gchi2_xz = new TH2D("h_gchi2_xz", Form("h_gchi2_xz; %s;%s", title_map[grid_x.f_name].c_str(), title_map[grid_z.f_name].c_str()), grid_x.f_N, grid_x.f_min, grid_x.f_max, grid_z.f_N, (grid_z.f_min)*2+1, (grid_z.f_max)*2+1);

		//minimize over two parameters
		TH1D* h_chi_delta = new TH1D("h_chi_delta", Form("h_chi_delta; %s;#Delta#chi^{2}",title_map[grid_z.f_name].c_str()), grid_z.f_N, (grid_z.f_min)*2+1, (grid_z.f_max)*2+1);

		//vectors that stores marginalized chi
		std::vector<double> mchi_xy, mchi_yz, mchi_xz, gchi_xy, gchi_yz, gchi_xz;

		for(int ix=1;ix <= grid_x.f_N; ix++){
		        for(int iy=1; iy <= grid_y.f_N; iy++) h_mchi2_xy->SetBinContent(ix, iy, DBL_MAX);
        		for(int iz=1; iz <= grid_z.f_N; iz++) h_mchi2_xz->SetBinContent(ix, iz, DBL_MAX);
   		}

   		for(int iz=1; iz <= grid_z.f_N; iz++){
        		for(int iy=1; iy <= grid_y.f_N; iy++)  h_mchi2_yz->SetBinContent(iy, iz, DBL_MAX);
        		h_chi_delta->SetBinContent(iz, DBL_MAX);
   		}


		for(int ix=0; ix < grid_x.f_N; ix++){
       		    for(int iy=0; iy < grid_y.f_N; iy++){
           		for(int iz=0 ; iz< grid_z.f_N; iz++){
                	    int ip = ix*grid_y.f_N*grid_z.f_N + iy*grid_z.f_N + iz; // index of grid point
		            std::vector<double> point = m_vec_grid[ip];


                            //marginalized minimum
                            //conditional operator, saver the smaller chi.
                	    if(vec_chi[ip]< h_mchi2_xy->GetBinContent(ix+1, iy+1)){
                        	 h_mchi2_xy->SetBinContent(ix+1, iy+1, vec_chi[ip]);
                         	//std::cout << "chi2 value: " << chi[ip] << std::endl;

		 	    }
                	    if(vec_chi[ip]< h_mchi2_xz->GetBinContent(ix+1, iz+1)) h_mchi2_xz->SetBinContent(ix+1, iz+1, vec_chi[ip]);
                	    if(vec_chi[ip]< h_mchi2_yz->GetBinContent(iy+1, iz+1)) h_mchi2_yz->SetBinContent(iy+1, iz+1, vec_chi[ip]);


               		    //global minimum
                	    if(point[2] == temp_best_point[2]){
				 h_gchi2_xy->Fill(point[0], point[1], vec_chi[ip]);
				 gchi_xy.push_back(vec_chi[ip]);}
                	    if(point[1] == temp_best_point[1]){
				 h_gchi2_xz->Fill(point[0], point[2]*2+1, vec_chi[ip]);
				 gchi_xz.push_back(vec_chi[ip]);}
               		    if(point[0] == temp_best_point[0]){
				 h_gchi2_yz->Fill(point[1], point[2]*2+1, vec_chi[ip]);
				 gchi_yz.push_back(vec_chi[ip]);}

                 	    //marginalize two parameters
                	    if(vec_chi[ip] < h_chi_delta->GetBinContent(iz+1)) h_chi_delta->SetBinContent(iz+1, vec_chi[ip]);
           		}
        	    }
   		}

		for(int iy=1;iy <= grid_y.f_N; iy++){
		        for(int ix=1; ix <= grid_x.f_N; ix++) mchi_xy.push_back(h_mchi2_xy->GetBinContent(ix, iy));
   		}

   		for(int iz=1; iz <= grid_z.f_N; iz++){
        		for(int iy=1; iy <= grid_y.f_N; iy++)  mchi_yz.push_back(h_mchi2_yz->GetBinContent(iy, iz));
        		for(int ix=1; ix <= grid_x.f_N; ix++)  mchi_xz.push_back(h_mchi2_xz->GetBinContent(ix, iz));
   		}
		
		h_mchi2_xy->Write(); h_mchi2_xz->Write(); h_mchi2_yz->Write();
		h_gchi2_xy->Write(); h_gchi2_xz->Write(); h_gchi2_yz->Write();
		h_chi_delta->Write();
	
		TCanvas c("c_chi_delta", "c_chi_delta");
		h_chi_delta->Draw("hist");
		TLine line(h_chi_delta->GetXaxis()->GetXmin(), 1.0, h_chi_delta->GetXaxis()->GetXmax(), 1.0);
		TLine line90(h_chi_delta->GetXaxis()->GetXmin(), 2.71, h_chi_delta->GetXaxis()->GetXmax(), 2.71);
		line.SetLineColor(8);
		line90.SetLineColor(42);
		line.Draw("same");
		line90.Draw("same");
		c.Update();
		c.SaveAs((tag+"_pretty_margin_chi.pdf").c_str(), "pdf");
		c.Write();



		if(grid_x.f_N >=5 && grid_y.f_N>=5){
			//draw 1,2,3 sigma contours for h_mchi2_zy
			TH2D h_mchi2_xy_inter =this->Do2DInterpolation(m_interpolation_number,grid_x.f_points, grid_y.f_points, mchi_xy, tag+"_XY"); 
			std::vector<TGraph*> g_mchi_contour = this->FindContour(h_mchi2_xy_inter, 2, tag+"_XY");

			h_mchi2_xy->GetXaxis()->SetRangeUser(h_mchi2_xy_inter.GetXaxis()->GetXmin(), h_mchi2_xy_inter.GetXaxis()->GetXmax());
			h_mchi2_xy->GetYaxis()->SetRangeUser(h_mchi2_xy_inter.GetYaxis()->GetXmin(), h_mchi2_xy_inter.GetYaxis()->GetXmax());
			std::cout << "x max: " << h_mchi2_xy->GetXaxis()->GetXmax() << " "<< h_mchi2_xy_inter.GetXaxis()->GetXmax() << std::endl;
			this->DrawContour(h_mchi2_xy, g_mchi_contour, std::vector<double>{});
			//TCanvas* c_2Dplot =(TCanvas*)this->DrawContour(h_mchi2_xy, g_mchi_contour, std::vector<double>{}).Clone();
			//fin->cd();
			//c_2Dplot->Write();
		}
	   } //end of the m_bool_poly_grid=false loop
	}//end of 3 dimension case
    fin->Close();
    } //end of map check

    return 0;
}


//2D interpolation
//x, y vector should only be in increasing order
//NOTE::elements of the vector should first change with vector 'x', then change with vector 'y'!! Order is important here!
TH2D SBNsinglephoton::Do2DInterpolation(int inter_number, std::vector<double>& x, std::vector<double>& y, std::vector<double>& value, std::string intag){
	double x_min = x[0];
	double x_max = x[x.size()-1];
	double y_min = y[0];
	double y_max = y[y.size() -1];
	//interpolation step size
	double x_step = fabs(x_max - x_min)/double(inter_number -1);
	double y_step = fabs(y_max - y_min)/double(inter_number -1);

	std::cout << "SBNsinglephoton::Do2DInterpolation\t||\tCreate interpolated 2D plot with tag: "<< intag << ", interpolation number " << inter_number<< std::endl;
	//delete gROOT->FindObject("h_inter");
	TH2D h_inter(("h_inter_"+intag).c_str(), ("h_inter_"+intag).c_str(), inter_number, x_min, x_max, inter_number, y_min, y_max);
	
	const gsl_interp2d_type *T = gsl_interp2d_bicubic;  //bicubic interpolation
	gsl_spline2d *spline = gsl_spline2d_alloc(T, x.size(), y.size());  
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
        gsl_interp_accel *yacc = gsl_interp_accel_alloc();

        /* initialize interpolation */
	gsl_spline2d_init(spline,  &x[0], &y[0], &value[0] , x.size(), y.size());
	for (int i = 0; i < inter_number; i++){

              double xi;

              if(i == (inter_number -1)) xi = x_min + (i-0.5)*x_step;  //to fill the last bin
              else xi = x_min + i*x_step;
	    
              for (int j = 0; j < inter_number; j++){

                  double yj;
                  if(j == (inter_number -1)) yj = y_min + (j-0.5)*y_step;
                  else yj = y_min + j*y_step;

		  double zij = gsl_spline2d_eval(spline, xi, yj, xacc, yacc);  
		  h_inter.Fill(xi, yj, zij);
	      }
	}

	//freee up the pointers
	gsl_spline2d_free(spline);
        gsl_interp_accel_free(xacc);
        gsl_interp_accel_free(yacc);
	

	return h_inter;	
}


//return 1/2/3/4/5 sigma contours as vector of TGraph
std::vector<TGraph*> SBNsinglephoton::FindContour(TH2D hin, int n, std::string intag){
   if(is_verbose) std::cout<< "SBNsinglephoton::FindContour\t||\tTrying to find " << n << "contours" << std::endl;
   
   std::vector<double> full_contour{2.30, 6.18,11.83, 19.35, 28.23}; // chi2 value for 2 dof with 1, 2, 3, 4, 5 sigma confidence level
   std::vector<std::string> full_ctour_string{"1sigma","2sigma", "3sigma", "4sigma", "5sigma"};
   std::vector<std::string> full_name_string{"1#sigma", "2#sigma", "3#sigma", "4#sigma", "5#sigma"};
   std::vector<int> full_color{kGreen-7, kCyan-7, kMagenta-7, kRed-7, kBlue-7};//color of the contour: kBlue, kMagenta, kRed
   std::vector<int> full_style{kSolid, kDashed, kDotted, kDashDotted, 5};
   std::vector<int> full_marker_style{24,25,26,27,28};

   std::vector<double> contour(full_contour.begin(), full_contour.begin()+n); // chi2 value for 2 dof with 1 sigma, 90%, 2 sigma and 3 sigma confidence level
   std::vector<std::string> CL_string(full_ctour_string.begin(), full_ctour_string.begin()+n);
   std::vector<std::string> Name_string(full_name_string.begin(), full_name_string.begin()+n);
   std::vector<int> color_vec(full_color.begin(), full_color.begin()+n); 
   std::vector<int> style_vec(full_style.begin(), full_style.begin()+n);
   std::vector<int> mstyle_vec(full_marker_style.begin(), full_marker_style.begin()+n); 

   //draw contour
   TCanvas c((intag+"c").c_str(), (intag+"c").c_str());
   c.cd();
   hin.SetContour((int)contour.size(), &contour[0]);
   hin.Draw("CONT Z LIST");//"LIST" generates a list of TGraph for each contour
   
   std::vector<TGraph*> out_graph;
   out_graph.resize(n);
   //TGraph vg_temp();
   //std::vector<TGraph> out_graph(n,vg_temp);	
   //std::cout << "check " << __LINE__ << std::endl;

   //grab contour object
   c.Update();
   std::cout << "check: " <<__LINE__ << std::endl;
   TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
   //loop over contours
   for(int i=0; i<conts->GetSize(); i++){
	//create a new TList pointer inside a loop, it will get deleted everytime during the loop --> no memory leak here.
        TList* con_list = (TList*)conts->At(i);   //a list of TGraph for i'th contour.
        double x_min = DBL_MAX;
        double x_max = DBL_MIN;
        double y_min = DBL_MAX;
        double y_max = DBL_MIN;

        out_graph[i] = new TGraph();   // one TGraph for one contour
        out_graph[i]->SetName(Form("graph_%s", CL_string[i].c_str()));
        out_graph[i]->SetTitle(Form("%s", Name_string[i].c_str()));
        //out_graph[i]->SetLineColor(color_vec[i]);
	out_graph[i]->SetLineStyle(style_vec[i]);
        out_graph[i]->SetLineWidth(2);
        out_graph[i]->SetMarkerColor(color_vec[i]);
	out_graph[i]->SetMarkerStyle(mstyle_vec[i]);
	out_graph[i]->SetMarkerSize(0.1);
        TGraph* c_graph = (TGraph*)con_list->First();  //grab the TGraph


        for(int j=0; j< con_list->GetSize() ; j++){
                TGraph* temp_graph= (TGraph*)c_graph->Clone();
		double x,y;
                for(int k =0; k< temp_graph->GetN(); k++){
                        temp_graph->GetPoint(k, x, y);
                        if(x < x_min) x_min = x;
                        if(x > x_max) x_max = x;

                        if(y < y_min) y_min =y;
                        if(y > y_max) y_max = y;

                        out_graph[i]->SetPoint(out_graph[i]->GetN(), x, y);
                }
		c_graph=(TGraph*)con_list->After(c_graph);
        }


	std::cout << "Contour " << CL_string[i] << ": " << contour[i] << std::endl;
        std::cout << "range for x : " << std::setprecision(3) << x_min << "~" << x_max << std::endl;
        std::cout << "range for y : " << std::setprecision(3) << y_min << "~" << y_max << std::endl;
   }

   //conts->SetOwner(kTRUE);
   //conts->Clear();
   return out_graph;
}



void SBNsinglephoton::DrawContour(TH2D* h, std::vector<TGraph*>& v_graph, std::vector<double> bf_point){
	this->DrawContour(h, v_graph, tag, bf_point);
	return;
}

void SBNsinglephoton::DrawContour(TH2D* h, std::vector<TGraph*>& v_graph, std::string intag,  std::vector<double> bf_point){

   if(is_verbose) std::cout << "SBNsinglephoton::DrawContour||\tCreating contour pdf: "<<intag<<"_pretty_chi_contours.pdf"<< std::endl;
   delete gROOT->FindObject((intag+"_pretty_contour").c_str());
   TCanvas c_canvas((intag+"_pretty_contour").c_str(), (intag+"_pretty_contour").c_str());
   gStyle->SetOptStat(0);
   //TLatex txt(3.8, 1.0, "#splitline{MicroBooNE}{Preliminary}");
   TLegend leg(0.68, 0.7, 0.9,0.9);
   leg.SetFillStyle(0); 
   leg.SetBorderSize(0);
   h->Draw("colz");
   TMarker* marker=nullptr;
   if(bf_point.size() == 2){
	marker=new TMarker(bf_point.at(1), bf_point.at(0), 29);
	marker->SetMarkerColor(kWhite);
	marker->SetMarkerSize(2);
	marker->Draw();
	leg.AddEntry(marker, "Best-fit Point", "P");
   }

   c_canvas.cd();
   for(auto& g:v_graph){
        g->SetLineColor(kWhite);
        g->SetMarkerColor(kWhite);
        g->Draw("same l");
        //g->Draw("same C");
        leg.AddEntry(g, g->GetTitle(), "L");
  }
  leg.SetTextSize(0.043);
  leg.Draw();
  c_canvas.Update();
  c_canvas.SaveAs((intag+"_pretty_chi_contours.pdf").c_str(),"pdf");
  delete marker;
  return;
}


int SBNsinglephoton::PrintOutFitInfo(std::map<int, std::vector<double>>& inmap, std::string intag, bool print_all){
	std::map<int, std::vector<double>> map_chi = inmap;
        std::map<int, std::vector<double>>::iterator itmap = map_chi.begin();

        int best_point = itmap->first;
        std::vector<double> vec_chi = itmap->second;
        double chi_min = *std::min_element(vec_chi.begin(), vec_chi.end());

	//best-fit point index at two grid
	int m_poly_grid_index = best_point/m_num_total_gridpoints;
	int m_grid_index = best_point%m_num_total_gridpoints;

	//print out info
	std::vector<double> m_point = m_vec_grid[m_grid_index];
	std::cout << intag<<": Best chi value is " << chi_min << " at flat point";
	for(int i=0; i<m_point.size(); i++)
	   std::cout << ", " << m_grid.f_dimensions[i].f_name << ": "<< m_point[i] ;
	if(m_bool_poly_grid){
	   std::vector<double> m_poly_point = m_vec_poly_grid[m_poly_grid_index];
	   
	   std::cout << ", and energy/momentum dependent shift for NCpi0 non-coherent component with";
	   for(int i=0 ; i< m_poly_point.size();i++)
		std::cout << ", "<< i+1 <<"nd order factor "<< m_poly_point[i];
	}
	std::cout << ". "<< std::endl;

	if(print_all){
	   std::cout << intag<< "========================Below are detailed chi2 and coordinate values==========================" << std::endl;
	   for(int i=0;i<vec_chi.size(); i++){
		int temp_poly_grid_index = i/m_num_total_gridpoints;
		int temp_grid_index = i%m_num_total_gridpoints;
		std::cout << "Point " << i << " Chi:" << vec_chi[i] << " Coordinate: ";
		std::vector<double> temp_point = m_vec_grid[temp_grid_index];
		for(int j=0 ;j<temp_point.size();j++) std::cout << temp_point[j] << ", ";
		if(m_bool_poly_grid){
		    std::vector<double> temp_poly_point = m_vec_poly_grid[temp_poly_grid_index];
		    std::cout << "PolyFactor: ";
		    for(int j=0; j<temp_poly_point.size();j++) std::cout << j+1<<"nd order: "<< temp_poly_point[j]<< ", ";
		}
		std::cout<< std::endl;
	   }
	   std::cout <<intag << "=========================End of detailed chi info=============================================="<< std::endl;
	}
	return 0;
}

int SBNsinglephoton::WriteOutInfo(std::map<int, std::vector<double>>& inmap){
	std::map<int, std::vector<double>>::iterator itmap = inmap.begin();

	TFile* fout = new TFile(Form("%s_fit_output.root", tag.c_str()), "UPDATE");
	TVectorD v(1);
	v[0] = itmap->first;
	std::vector<double> vchi=itmap->second;

	fout->cd();
	v.Write("bf_index");
	fout->WriteObject(&vchi, "vector_chi");
	fout->Close();
	return 0;
}


int SBNsinglephoton::SetFlatFullFracCovarianceMatrix(double flat){

	m_full_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
	for(int i=0;i<num_bins_total; i++) (*m_full_fractional_covariance_matrix)(i,i) = flat*flat;

	return 0;	
}

int SBNsinglephoton::SetStatOnly(){
	m_full_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
	m_full_but_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
	m_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);

	m_full_fractional_covariance_matrix->Zero();
	m_full_but_genie_fractional_covariance_matrix->Zero();
	m_genie_fractional_covariance_matrix->Zero();

	return 0;
}

int SBNsinglephoton::ModifyCV(double infactor){

	std::vector<double> temp_constrain_param;
	for(int i=0 ; i<m_grid.f_num_dimensions; i++){
	    if(m_grid.f_dimensions[i].f_has_constrain){
		temp_constrain_param.push_back(m_grid.f_dimensions[i].f_constrain_value);
	    }
	}

	if(m_bool_poly_grid){

           for(int i=0; i< m_poly_grid.f_num_dimensions; i++){
              if(m_poly_grid.f_dimensions[i].f_has_constrain){
		     temp_constrain_param.push_back(m_poly_grid.f_dimensions[i].f_constrain_value);
              }
	   }
	}
	
	return this->ModifyCV(infactor, temp_constrain_param);
}

int SBNsinglephoton::ModifyCV(double infactor, std::vector<double> param){


	if(!m_bool_cv_spectrum_loaded){
                this->LoadCV();
        }


	std::cout<< "SBNsinglephoton::ModifyCV\t|| Start modifying CV spectrum" <<std::endl;
	int index = 0;
	//apply flat normalization
	for(int i=0 ; i<m_grid.f_num_dimensions; i++){
            if(m_grid.f_dimensions[i].f_has_constrain  && (index < param.size()) ){
                m_cv_spectrum->Scale(m_grid.f_dimensions[i].f_name, param[index]);
		index++;
            }
        }

	//energy/momentum dependent scaling
	//now that we reset negative weight to 1, this needs to be modified, tho it's probably never gonna be used
	if(m_bool_poly_grid){
	   std::ostringstream prescale_tag;
	   std::vector<double> temp_scale_parameter;
           for(int i=0; i< m_poly_grid.f_num_dimensions; i++){
              if(m_poly_grid.f_dimensions[i].f_has_constrain && (index < param.size())  ){   
                     prescale_tag << "_" << std::fixed<< std::setprecision(3) << param[index];
		     temp_scale_parameter.push_back(param[index]);
		     index++;
              }
           }

	   if(temp_scale_parameter.size() != 0){
	       std::string temp_filename = tag+"_PreScaled"+prescale_tag.str()+".SBNspec.root";
	       //check if a file exits
	       if(gSystem->AccessPathName(temp_filename.c_str())){
		   std::cout << "SBNsinglephoton::ModifyCV\t|| Prescaled file doesn't exist, start generating it..." << std::endl;
		   if(!m_file_open) this->OpenFiles();
		   this->PreScaleSpectrum(xmlname, temp_scale_parameter);
	       }

	       std::cout <<"SBNsinglephoton::ModifyCV\t|| Adding pre-scaled spectrum to the CV." << std::endl;
	       SBNspec temp_prescale(temp_filename.c_str(), xmlname, false);
	       m_cv_spectrum->Add(&temp_prescale); 
	   }
	}

	if(m_bool_data_spectrum_loaded){
		 m_cv_spectrum->Scale("NCDeltaRadOverlayLEE", (infactor-1)*0.5);
		 m_cv_spectrum->CompareSBNspecs(m_data_spectrum, "ModifiedCV_Data");
        }

	m_bool_modify_cv = true;
	m_cv_delta_scaling = infactor;

	m_cv_spectrum->WriteOut("ModifiedCV");

	return 0;
}

SBNspec SBNsinglephoton::GeneratePointSpectra(int np){
	int m_poly_grid_index = np/m_num_total_gridpoints;
        int m_grid_index = np%m_num_total_gridpoints;

	//SBNspec spec_rgrid_part(tag+"_CV.SBNspec.root", xmlname, false);
	SBNspec spec_rgrid_part = *m_cv_spectrum;

//	spec_rgrid_part.Scale("NCPi0Coh", 2.0);
//        spec_rgrid_part.Scale("NCPi0NotCoh", 0.8);

	std::vector<double> temp_p_grid = m_vec_grid[m_grid_index];
	for(int i=0;i<m_grid.f_num_dimensions; i++){
	    if(m_grid.f_dimensions[i].f_name == "NCDeltaRadOverlayLEE" && temp_p_grid[i]<0 ){
		spec_rgrid_part.Scale("NCDeltaRadOverlayLEE", 0.0);
		spec_rgrid_part.Scale("NCDeltaRadOverlaySM", 1+temp_p_grid[i]*2.0);
	    }
	    else
		spec_rgrid_part.Scale(m_grid.f_dimensions[i].f_name, temp_p_grid[i]);
		//spec_rgrid_part.ScaleAll(temp_p_grid[i]);
	}

	if(m_bool_poly_grid){
	    double non_coh_factor;   
	    for(int i=0;i<m_grid.f_num_dimensions; i++){
		if(m_grid.f_dimensions[i].f_name == "NCPi0NotCoh"){
		    non_coh_factor=temp_p_grid[i];
		    spec_rgrid_part.Scale("NCPi0NotCoh", 0.0);
		}
	    }

	    std::vector<double> temp_p_polygrid=m_vec_poly_grid[m_poly_grid_index];

	    if(!m_file_open) this->OpenFiles();
	    if(!m_bool_cv_spectrum_generated) m_bool_cv_spectrum_generated=true;
	    SBNspec spec_polygrid_part(xmlname, -1, false);
	    this->ScaleSpectrum(&spec_polygrid_part, non_coh_factor, temp_p_polygrid);

	    spec_rgrid_part.Add(&spec_polygrid_part);
	}


	return spec_rgrid_part;
}


int SBNsinglephoton::PoissonFluctuation(SBNspec *inspec){
	std::cout << "SBNsinglephoton::PoissonFluctuation\t || Poisson fluctuation the spectrum!!" << std::endl;
	inspec->CollapseVector();

	inspec->ScalePoisson();
	inspec->CollapseVector();

	return 0;
}


double SBNsinglephoton::CalcChi(bool use_cnp){

	m_chi = new SBNchi(*m_cv_spectrum, m_full_fractional_covariance_matrix);
	if(!use_cnp){
	    return m_chi->CalcChi(m_data_spectrum);
	}else{
	    m_cv_spectrum->CollapseVector();
	    m_data_spectrum->CollapseVector();

	    TMatrixT<double> Mtemp = m_chi->CalcCovarianceMatrixCNP(*m_full_fractional_covariance_matrix, m_cv_spectrum->full_vector, m_cv_spectrum->collapsed_vector, m_data_spectrum->collapsed_vector);
	    TMatrixT<double> Invert_temp = m_chi->InvertMatrix(Mtemp);
	    return m_chi->CalcChi(Invert_temp, m_cv_spectrum->collapsed_vector, m_data_spectrum->collapsed_vector, true);
	}

}

void SBNsinglephoton::SetInterpolationNumber(int in){
	m_interpolation_number = in;
}
