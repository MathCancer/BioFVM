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

#include "../BioFVM.h" 

using namespace BioFVM;

int omp_num_threads = 8; // set number of threads for parallel computing
// set this to # of CPU cores x 2 (for hyperthreading)
double pi= 3.1415926535897932384626433832795;
double max_substrate_1=38.0;

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
					(*write_here)[0] = max_substrate_1;
				}	
				
	return ;
}

void supply_target_function( Microenvironment* m, int voxel_index , std::vector<double>* write_here )
{
	(*write_here)[0] = max_substrate_1; 	
	return ;
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

int write_report(int num_cells, double time)
{
	std::ofstream report ("scaling_test_results_numcells.txt", std::ofstream::app);
	report<<num_cells<<"\t"<<time<<std::endl;
	report.close();
	return 0;
}

int main( int argc, char* argv[] )
{
	double t = 0.0; 
	double t_output_interval = 1.0; 
	double t_next_output_time = 0;
	int next_output_index = 0;
	
	double cell_radius=5.0;
	double num_cells= strtod(argv[1],NULL);
	
	double spacing= pow(num_cells,0.333333)/800.0;  //1000-(40 micron strip + 60 micron empty space)*2 =800
	double dt = 0.01;
	double t_max = 2.0;
	double mesh_resolution=20;
	double cell_uptake_rate=10.0;
	double cell_secretion_rate=10.0;
	// openmp setup
	omp_set_num_threads(omp_num_threads);
		
	// create a microenvironment; 
	Microenvironment microenvironment;
	microenvironment.set_density(0, "substrate1" , "dimensionless" );
	
	std::vector<double> center(3);
	center[0] = 0; center[1] = 0; center[2] = 0; 
	
	double minX=-500, minY=-500, minZ=-500, maxX=500, maxY=500, maxZ=500;
	microenvironment.resize_space_uniform( minX,maxX,minY,maxY,minZ,maxZ, mesh_resolution );
	// register the diffusion solver 	
	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
	microenvironment.diffusion_coefficients[0] = 100000; // microns^2 / min 
	microenvironment.decay_rates[0] = 0.1;
	microenvironment.display_information( std::cout );

	microenvironment.bulk_supply_rate_function =  supply_function;
	microenvironment.bulk_supply_target_densities_function = supply_target_function;
	
	#pragma omp parallel for 
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		microenvironment.density_vector(i)[0]= max_substrate_1; 
	}
	
	std::vector<double> tempPoint(3,0.0);
	int counter=0;
	for(double x=-400;x<400;x+=spacing)
		for(double y=-400;y<400;y+=spacing)
			for(double z=-400;z<400;z+=spacing)
			{
				if(counter>=num_cells)
					break;
				tempPoint[0]=x;
				tempPoint[1]=y;
				tempPoint[2]= z;
								
				Basic_Agent * temp_point_sink = create_basic_agent();
				temp_point_sink->register_microenvironment(&microenvironment);
				temp_point_sink->assign_position(tempPoint);
				temp_point_sink->volume= (4.0/3.0)*pi*pow(cell_radius,3.0);
				(*temp_point_sink->uptake_rates)[0]=cell_uptake_rate;
				temp_point_sink->set_internal_uptake_constants(dt); 
				counter++;
			}
			
	// display information 
	microenvironment.display_information( std::cout );
	int center_voxel_index = microenvironment.nearest_voxel_index(center);
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
		microenvironment.simulate_cell_sources_and_sinks( dt );
		microenvironment.simulate_diffusion_decay( dt );
		t += dt; 
		output_index++; 
	}
	process_output(t_max, dt, mesh_resolution);
	write_report(num_cells, BioFVM::runtime_stopwatch_value());
	std::cout << "done!" << std::endl; 
	return 0; 
}