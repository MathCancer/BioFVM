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
double substrate1_supply_rate=10;
double substrate1_saturation_density=1;
double substrate1_uptake_rate=10;

void uptake_function( Microenvironment* microenvironment, int voxel_index , std::vector<double>* write_here )
{
	double spheroid_radius=100;
	std::vector<double> center(3);
	center[0] = microenvironment->mesh.bounding_box[0]+microenvironment->mesh.bounding_box[3];
	center[1] = microenvironment->mesh.bounding_box[1]+microenvironment->mesh.bounding_box[4];
	center[2] = microenvironment->mesh.bounding_box[2]+microenvironment->mesh.bounding_box[5];
	if(sqrt( norm_squared(microenvironment->voxels(voxel_index).center - center))<spheroid_radius)
		(*write_here)[0] = substrate1_uptake_rate; 	
	return ;
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
	
	#pragma omp parallel for 
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		microenvironment.density_vector(i)[0]= 1.0; 
	}
	
	// add dirichlet nodes 
	std::vector<double> dirichlet_zero( 1 , 1.0 );
	double min_x=microenvironment.mesh.bounding_box[0];
	double max_x=microenvironment.mesh.bounding_box[3];
	double min_y=microenvironment.mesh.bounding_box[1];
	double max_y=microenvironment.mesh.bounding_box[4];
	double min_z=microenvironment.mesh.bounding_box[2];
	double max_z=microenvironment.mesh.bounding_box[5];
	double strip_width=40;	

	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		if( abs(max_x-microenvironment.voxels(i).center[0]) < strip_width || abs(microenvironment.voxels(i).center[0]- min_x)< strip_width  
			|| abs(max_y-microenvironment.voxels(i).center[1]) < strip_width || abs(microenvironment.voxels(i).center[1]- min_y)< strip_width  
				|| abs(max_z-microenvironment.voxels(i).center[2]) < strip_width || abs(microenvironment.voxels(i).center[2]- min_z)< strip_width )
				{
					microenvironment.add_dirichlet_node( i , dirichlet_zero );
				}		
	}

	microenvironment.write_to_matlab( "initial_concentration.mat" );

	// register the diffusion solver 	
	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
	
	// register substrates properties 
	microenvironment.diffusion_coefficients[0] = 1000; // microns^2 / min 
	microenvironment.decay_rates[0] = 0.01;
				
	// display information 	
	microenvironment.display_information( std::cout );
	
	// register defined supply/uptake functions with the solver
	microenvironment.bulk_uptake_rate_function = uptake_function;
	
	double dt = 0.01; 	
	double t = 0.0; 
	double t_max=5;

	while( t < t_max )
	{
		microenvironment.simulate_bulk_sources_and_sinks( dt );
		microenvironment.simulate_diffusion_decay( dt );
		t += dt; 
	}
	microenvironment.write_to_matlab( "final.mat" );
	std::cout<<"done!"<<std::endl;
	return 0; 
}