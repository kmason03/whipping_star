#Completely stolen.
CPPFLAGS=-I $(BOOST_INC) \
         -I $(CANVAS_INC) \
         -I $(CETLIB_INC) \
         -I $(FHICLCPP_INC) \
         -I $(LARCOREOBJ_INC) \
         -I $(LARDATAOBJ_INC) \
	 -I $(LARSIM_INC) \
         -I $(NUSIMDATA_INC) \
	 -I $(UBOONECODE_INC) \
         -I $(ROOT_INC)

# uboonecode v08_00_00_48
UBOONECODE_LIB="/cvmfs/uboone.opensciencegrid.org/products/uboonecode/v08_00_00_48/slf7.x86_64.e17.prof/lib"

CXXFLAGS=-std=c++17 -Wall
CXX=g++
LDFLAGS=$$(root-config --libs) \
        -L $(CANVAS_LIB) -l canvas \
        -L $(CETLIB_LIB) -l cetlib \
        -L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
        -L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
        -L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase \
        -L $(LARSIM_LIB) -l larsim_EventWeight_Base


UNAME := $(shell uname -s)

all: form_input_to_sbnfit_run1 form_input_to_sbnfit_run2 form_input_to_sbnfit_run3 form_input_to_sbnfit_forJosh1 form_input_to_sbnfit_forJosh2 form_input_to_sbnfit_forJosh3

form_input_to_sbnfit_run1: form_input_to_sbnfit_forAppearence_run1.cpp
	@echo Building $(EXEC)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

form_input_to_sbnfit_run2: form_input_to_sbnfit_forAppearence_run2.cpp
	@echo Building $(EXEC)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

form_input_to_sbnfit_run3: form_input_to_sbnfit_forAppearence_run3.cpp
	@echo Building $(EXEC)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

form_input_to_sbnfit_forJosh1: form_input_to_sbnfit_forJosh1.cpp
	@echo Building $(EXEC)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

form_input_to_sbnfit_forJosh2: form_input_to_sbnfit_forJosh2.cpp
	@echo Building $(EXEC)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

form_input_to_sbnfit_forJosh3: form_input_to_sbnfit_forJosh3.cpp
	@echo Building $(EXEC)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
clean: 
	rm -f form_input_to_sbnfit
