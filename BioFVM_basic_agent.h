/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.0.4) [1]        #
#                                                                           #
# [1] A. Ghaffarizaeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient  #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics, 2015. DOI: 10.1093/bioinformatics/btv730               #
#############################################################################
#                                                                           #
# Copyright 2015 Paul Macklin and the BioFVM Project                        #
#                                                                           #
# Licensed under the Apache License, Version 2.0 (the "License");           #
# you may not use this file except in compliance with the License.          #
# You may obtain a copy of the License at                                   #
#                                                                           #
#    http://www.apache.org/licenses/LICENSE-2.0                             #
#                                                                           #
# Unless required by applicable law or agreed to in writing, software       #
# distributed under the License is distributed on an "AS IS" BASIS,         #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #
# See the License for the specific language governing permissions and       #
# limitations under the License.                                            #
#############################################################################
*/

#ifndef __BioFVM_basic_agent_h__
#define __BioFVM_basic_agent_h__

#include <vector>
#include "BioFVM_microenvironment.h"
#include "BioFVM_matlab.h"
#include "BioFVM_vector.h"

namespace BioFVM{


class Basic_Agent
{
 private:
	Microenvironment* microenvironment; 
	int selected_microenvironment; 
	std::vector<double> cell_source_sink_solver_temp1;
	std::vector<double> cell_source_sink_solver_temp2;
	int current_microenvironment_voxel_index;
	
 protected:
	std::vector<double> previous_velocity; 
	int current_voxel_index;	
 public:
	std::vector<double> * secretion_rates; 
	std::vector<double> * saturation_densities; 
	std::vector<double> * uptake_rates;  
	
	void set_internal_uptake_constants( double dt ); // any time you update the cell volume or rates, should call this function. 

	void register_microenvironment( Microenvironment* );
	Microenvironment* get_microenvironment( void ); 

	int ID; 
	int index; 
	int type;
	double volume;
	
	bool assign_position(double x, double y, double z);
	bool assign_position(std::vector<double> new_position);
	
	std::vector<double> position;  
	std::vector<double> velocity; 
	void update_position( double dt );
	
	Basic_Agent(); 
	~Basic_Agent(); 
	// simulate secretion and uptake at the nearest voxel at the indicated microenvironment.
	// if no microenvironment indicated, use the currently selected microenvironment. 
	void simulate_secretion_and_uptake( Microenvironment* M, double dt ); 

	// find the nearest voxel index for the indicated microenvironment 
	int nearest_voxel_index( int microenvironment_index ); 
	int nearest_voxel_index( void ); 
	// directly access the substrate vector at the nearest voxel at the indicated microenvironment 
	std::vector<double>& nearest_density_vector( int microenvironment_index ); 
	std::vector<double>& nearest_density_vector( void );
};

extern std::vector<Basic_Agent*> all_basic_agents; 

Basic_Agent* create_basic_agent( void );
void delete_basic_agent( int ); 
void delete_basic_agent( Basic_Agent* ); 
void save_all_basic_agents_to_matlab( std::string filename ); 

};

#endif

