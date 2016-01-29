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

int main( int argc, char* argv[] )
{
	double t = 0.0; 
	double t_max = strtod(argv[1],NULL);
	double mesh_resolution=strtod(argv[2],NULL);
	t=t_max;
	double dx, dy, dz; 
	
	// openmp setup
	omp_set_num_threads(omp_num_threads);
	
	// create a microenvironment; 
	Microenvironment microenvironment;
	microenvironment.name="substrate scale";

	microenvironment.set_density(0, "substrate1" , "dimensionless" );
	microenvironment.spatial_units = "microns";
	microenvironment.mesh.units = "microns";
	microenvironment.time_units = "minutes";
	
	std::vector<double> center(3);

	double minX=-500, minY=-mesh_resolution/2, minZ=-mesh_resolution/2, maxX=500, maxY=mesh_resolution/2, maxZ=mesh_resolution/2; 
	microenvironment.resize_space_uniform( minX,maxX,minY,maxY,minZ,maxZ, mesh_resolution );
	microenvironment.display_information( std::cout );
	center[0] = 0; center[1] = 0; center[2] = 0;
	int center_voxel_index = microenvironment.nearest_voxel_index(center);
	double t0=0.05;
	
	try {

	double L0=500.0, D=100000;
	// #pragma omp parallel for 
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		double new_value;
		double X = sqrt(norm_squared(microenvironment.voxels(i).center ));
		new_value=1.0+ cos((pi/L0) *X)* exp(-(D* pi*pi/ (L0*L0)) *t);
		microenvironment.density_vector(i)[0]=new_value;
	}
	microenvironment.write_to_matlab( "final.mat" );
	
	std::cout << "done!" << std::endl; 
	}
	catch( const std::exception& e ) { 
		std::cout << e.what(); 
	}
	return 0; 
}