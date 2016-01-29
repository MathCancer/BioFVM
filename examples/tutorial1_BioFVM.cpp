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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <random>

#include "../BioFVM.h" 

using namespace BioFVM;

int omp_num_threads = 8; // set number of threads for parallel computing
// set this to # of CPU cores x 2 (for hyperthreading)
double pi= 3.1415926535897932384626433832795;

double UniformRandom()
{	
	return ((double) rand() / (RAND_MAX));
}

int main( int argc, char* argv[] )
{	
	// openmp setup
	omp_set_num_threads(omp_num_threads);
	
	// create a microenvironment; 
	Microenvironment microenvironment;
	microenvironment.name="substrate scale";

	microenvironment.set_density(0, "substrate1" , "dimensionless" );
	microenvironment.spatial_units = "microns";
	microenvironment.mesh.units = "microns";
	microenvironment.time_units = "minutes";
	
	
	double minX=0, minY=0, minZ=0, maxX=1000, maxY=1000, maxZ=1000, mesh_resolution=10;
	microenvironment.resize_space_uniform( minX,maxX,minY,maxY,minZ,maxZ, mesh_resolution );
	microenvironment.display_information( std::cout );

	
	std::vector<double> center(3);
	center[0] = (microenvironment.mesh.bounding_box[0]+microenvironment.mesh.bounding_box[3])/2;
	center[1] = (microenvironment.mesh.bounding_box[1]+microenvironment.mesh.bounding_box[4])/2;
	center[2] = (microenvironment.mesh.bounding_box[2]+microenvironment.mesh.bounding_box[5])/2;
	double stddev_squared = -100.0 * 100.0; 	
	std::vector<double> one( microenvironment.density_vector(0).size() , 1.0 ); 
	#pragma omp parallel for 
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		std::vector<double> displacement = microenvironment.voxels(i).center - center; 
		double distance_squared = norm_squared( displacement );
		double coeff = distance_squared;
		coeff /=  stddev_squared;
		microenvironment.density_vector(i)[0]= exp( coeff ); 
	}
	microenvironment.write_to_matlab( "initial_concentration.mat" ); 

	
	

	// register the diffusion solver 	
	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
	
	// register substrates properties 
	microenvironment.diffusion_coefficients[0] = 1000; // microns^2 / min 
	microenvironment.decay_rates[0] = 0.01;
				
	// display information 	
	microenvironment.display_information( std::cout );
	
	double dt = 0.01; 
	double cell_radius=5;
	for(int i=0; i< 500;i++)
	{
		std::vector<double> tempPoint(3,0.0);
		for( int j=0; j < 3 ; j++ )
		{ tempPoint[j] = UniformRandom()*1000; }		
		
		Basic_Agent * temp_point_source = create_basic_agent();
		temp_point_source->register_microenvironment(&microenvironment);
		temp_point_source->assign_position(tempPoint);
		temp_point_source->volume= (4.0/3.0)*pi*pow(cell_radius,3.0);
		(*temp_point_source->secretion_rates)[0]=10;
		(*temp_point_source->saturation_densities)[0]=1;
		temp_point_source->set_internal_uptake_constants(dt); 
		

		for( int j=0; j < 3 ; j++ )
		{ tempPoint[j] = UniformRandom()*1000; }		
		Basic_Agent * temp_point_sink = create_basic_agent();
		temp_point_sink->register_microenvironment(&microenvironment);
		temp_point_sink->assign_position(tempPoint);
		temp_point_sink->volume= (4.0/3.0)*pi*pow(cell_radius,3.0);
		(*temp_point_sink->uptake_rates)[0]=0.8;
		temp_point_sink->set_internal_uptake_constants(dt); 
	}
	
	
	double t = 0.0; 
	double t_max=5;

	while( t < t_max )
	{
		microenvironment.simulate_cell_sources_and_sinks( dt );
		microenvironment.simulate_diffusion_decay( dt );
		t += dt; 
	}
	microenvironment.write_to_matlab( "final.mat" );
	std::cout<<"done!"<<std::endl;
	return 0; 
}