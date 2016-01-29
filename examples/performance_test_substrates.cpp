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

//setup Mersenne Twister random generator 
std::random_device rd;
std::mt19937 gen(rd());

std::vector<double> substrates_sources;
int num_substrates;

void supply_function( Microenvironment* microenvironment, int voxel_index , std::vector<double>* write_here )
{
	double max_xyz=500;
	double min_xyz=-500;
	double strip_width=40;
	(*write_here)[0] = 0; 

	if( abs(max_xyz-microenvironment->voxels(voxel_index).center[0]) < strip_width || abs(microenvironment->voxels(voxel_index).center[0]- min_xyz)< strip_width  
			|| abs(max_xyz-microenvironment->voxels(voxel_index).center[1]) < strip_width || abs(microenvironment->voxels(voxel_index).center[1]- min_xyz)< strip_width  
				|| abs(max_xyz-microenvironment->voxels(voxel_index).center[2]) < strip_width || abs(microenvironment->voxels(voxel_index).center[2]- min_xyz)< strip_width )
				{
					// std::cout<<"test"<<std::endl;
					for(int i=0;i< num_substrates;i++)
						(*write_here)[i] = substrates_sources[i];
				}	
				
	return ;
}

void supply_target_function( Microenvironment* m, int voxel_index , std::vector<double>* write_here )
{
	for(int i=0;i< num_substrates;i++)
		(*write_here)[i] = substrates_sources[i];	
	return ;
}

int write_report(double time)
{
	std::ofstream report ("scaling_test_results_substrates.txt", std::ofstream::app);
	report<<num_substrates<<"\t"<<time<<std::endl;
	report.close();
	return 0;
}

double get_rand(double lb, double ub)
{
	double rnd= std::generate_canonical<double, 10>(gen);;
	double range=ub-lb;
	return range*rnd+lb;
}

void process_output(double t, double dt, double mesh_resolution)
{

	std::cout << "current simulated time: " << t   << " minutes " << std::endl; 
	std::cout << "interval wall time: ";
	BioFVM::TOC();
	BioFVM::display_stopwatch_value( std::cout , BioFVM::stopwatch_value() ); 
	std::cout << std::endl; 
	std::cout << "total wall time: "; 
	BioFVM::RUNTIME_TOC();
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 
	std::cout << std::endl;
	
	std::cout << "time: "<<t<<std::endl;
	BioFVM::TIC();
}

int main( int argc, char* argv[] )
{
	double t = 0.0; 
	double t_output_interval = 1.0; 
	double t_next_output_time = 0;
	int next_output_index = 0;
	
	double dt = 0.01;
	double t_max = 2.0; 
	double mesh_resolution=10.0;
	
	num_substrates=atoi(argv[1]);
	
	substrates_sources.resize(num_substrates);
	// openmp setup
	omp_set_num_threads(omp_num_threads);
	
	// set random seed
	srand(2); 
	
	// create a microenvironment; 
	Microenvironment microenvironment;
	microenvironment.set_density(0, "substrate0" , "dimensionless" );
	microenvironment.diffusion_coefficients[0] = 100000; // microns^2 / min 
	microenvironment.decay_rates[0] = 0.1;
	
	 
	for(int i=1;i<num_substrates;i++)
	{
		std::string substrate_name; 
		substrate_name.resize( 25 , '\0' );
		sprintf( (char*) substrate_name.c_str() , "substrate%d" ,i ); 
		microenvironment.add_density(substrate_name , "dimensionless" );
		microenvironment.diffusion_coefficients[i] = get_rand(75000.0, 125000.0); // microns^2 / min 
		microenvironment.decay_rates[i] = get_rand(0.05, 0.15);
		substrates_sources.push_back(get_rand(30.0,60.0));
	}
	double minX=-500, minY=-500, minZ=-500, maxX=500, maxY=500, maxZ=500;//, mesh_resolution=10; 
	microenvironment.resize_space_uniform( minX,maxX,minY,maxY,minZ,maxZ, mesh_resolution );
	// register the diffusion solver 	
	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
	
	microenvironment.bulk_supply_rate_function =  supply_function;
	microenvironment.bulk_supply_target_densities_function = supply_target_function;
	
	#pragma omp parallel for 
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		for(int j=0;j< num_substrates;j++)
			microenvironment.density_vector(i)[j]= substrates_sources[j]; 
	}
	
	// display information 
	microenvironment.display_information( std::cout );
	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	int output_index=0;	
	
	std::cout<<"number of agents: "<< all_basic_agents.size() <<std::endl;

	while( t < t_max )
	{		
		if(  fabs( t - t_next_output_time ) < dt/10.0 )
		{
			process_output(t, dt, mesh_resolution);
			t_next_output_time += t_output_interval; 
			next_output_index++; 
		}
		// simulate microenvironment 
		microenvironment.simulate_bulk_sources_and_sinks( dt );
		microenvironment.simulate_diffusion_decay( dt );
		
		t += dt; 
		output_index++; 
	}
	process_output(t_max, dt, mesh_resolution);
	BioFVM::RUNTIME_TOC();
	write_report(BioFVM::runtime_stopwatch_value());
	std::cout << "done!" << std::endl; 
	return 0; 
}