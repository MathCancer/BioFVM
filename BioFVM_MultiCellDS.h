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

#ifndef __BioFVM_MultiCellDS_h__
#define __BioFVM_MultiCellDS_h__

#include "pugixml.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include <fstream> 

namespace BioFVM{
extern std::string MultiCellDS_version_string; 
extern std::string MultiCellDS_experimental_snapshot_type_string; 
extern int MultiCellDS_experimental_snapshot_code; 

extern std::string MultiCellDS_simulation_snapshot_type_string; 
extern int MultiCellDS_simulation_snapshot_code; 

extern std::string MultiCellDS_digital_cell_line_type_string; 
extern int MultiCellDS_digital_cell_line_code; 

class Microenvironment; 

extern pugi::xml_document biofvm_doc; 

class MultiCellDS_Metadata
{
	private:
	public:
		std::string MultiCellDS_type; 
	
		std::string user; 
		std::string program; 
		std::string program_version; 
		
		std::vector<double> bounding_box; 
		std::string spatial_units; 
		std::string time_units;
		std::string runtime_time_units; 
		double current_time; 
		double current_runtime; 
	
		MultiCellDS_Metadata();	
		void add_to_open_MultiCellDS_xml( std::ostream& os , std::string base_tabbing ); 
};

extern MultiCellDS_Metadata default_metadata; 

// functions to read multiscale_microenvironment from MultiCellDS file (requires pugixml)
void read_microenvironment_from_MultiCellDS_xml( Microenvironment& M_destination , std::string filename );
void read_microenvironment_from_MultiCellDS_xml( Microenvironment& M_destination , pugi::xml_document& xml_dom ); 

void set_save_biovm_mesh_as_matlab( bool newvalue ); 
void set_save_biovm_data_as_matlab( bool newvalue ); 

void save_microenvironment_to_MultiCellDS_xml_pugi( std::string filename_base , Microenvironment& M );
};

#endif 
