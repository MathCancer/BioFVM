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

#include "BioFVM.h"

#include <omp.h>

using namespace BioFVM; 

void supply_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
	// use this syntax to access the jth substrate write_here
	// (*write_here)[j]
	// use this syntax to access the jth substrate in voxel voxel_index of microenvironment: 
	// microenvironment->density_vector(voxel_index)[j]

	return; 
}

void supply_target_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
	// use this syntax to access the jth substrate write_here
	// (*write_here)[j]
	// use this syntax to access the jth substrate in voxel voxel_index of microenvironment: 
	// microenvironment->density_vector(voxel_index)[j]

	return; 
}

void uptake_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
	// use this syntax to access the jth substrate write_here
	// (*write_here)[j]
	// use this syntax to access the jth substrate in voxel voxel_index of microenvironment: 
	// microenvironment->density_vector(voxel_index)[j]

	return; 
}

int main( int argc, char* argv[] )
{
	omp_set_num_threads( 8 );
	
	std::cout << "Starting program ... " << std::endl;
	
	// create a microenvironment, and set units 
		
	Microenvironment M; 
	M.name = "microenvironment"; 
	M.time_units = "min"; 
	M.spatial_units = "micron"; 
	M.mesh.units = M.spatial_units;
	
	// set up and add all the densities you plan 

	M.set_density( 0 , "unnamed 1" , "dimensionless" ); 
	
	// here's how you add a new substrate 
	M.add_density( "unnamed 2" , "dimensionless" ); 
	
	
	// set the properties of the diffusing substrates 
	
	M.diffusion_coefficients[0] = 1e5;   
	M.decay_rates[0] = 10; // 100 micron length scale 
	M.diffusion_coefficients[1] = 1e5;   
	M.decay_rates[1] = 10; // 100 micron length scale 
	
	// set the mesh size 
	
	double dx = 20; // 
	M.resize_space( 0.0 , 1000.0 , 0, 1000.0 , 0.0 , 1000.0 , dx, dx, dx );  
	
	// display summary information 
	
	M.display_information( std::cout ); 

	// set initial conditions 
	
	// use this syntax to create a zero vector of length 3
	// std::vector<double> zero(3,0.0); 
	
	// use this syntax for a parallelized loop over all the 
	// voxels in your mesh: 	
	#pragma omp parallel for 
	for( int i=0 ; i < M.number_of_voxels() ; i++ )
	{
		// use this syntax to access the coordinates (as a vector) of 
		// the ith voxel; 
		// M.mesh.voxels[i].center 
		
		// use this access the jth substrate at the ith voxel
		// M.density_vector(i)[j]
		
	}

	// save the initial profile 
	
	M.write_to_matlab( "initial.mat" );
	
	// set up the diffusion solver, sources and sinks 
	
	M.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D;
	
	M.bulk_supply_rate_function = supply_function;
	M.bulk_supply_target_densities_function = supply_target_function;
	M.bulk_uptake_rate_function = uptake_function;
	
	double t     = 0.0; 
	double t_max = 100.0;
	double dt    = 0.1; 
	
	double output_interval  = 10.0;  // how often you save data 
	double next_output_time = t;     // next time you save data 
	
	while( t < t_max )
	{
		// if it's time, save the simulation 
		if( fabs( t - next_output_time ) < dt/2.0 )
		{
			std::cout << "simulation time: " << t << " " << M.time_units << " (" << t_max << " " << M.time_units << " max)" << std::endl; 
			
			char* filename; 
			filename = new char [1024];
			sprintf( filename, "output_%6f.mat" , next_output_time ); 
			M.write_to_matlab( filename ); 
			delete [] filename; 
			next_output_time += output_interval; 
		}
		
		M.simulate_bulk_sources_and_sinks( dt );
		M.simulate_diffusion_decay( dt );
		M.simulate_cell_sources_and_sinks( dt ); 
		
		t += dt;
	}
	
	M.write_to_matlab( "final.mat");
	
	std::cout << "Done!" << std::endl; 
	
	return 0; 
}
