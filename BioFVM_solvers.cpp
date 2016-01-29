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

#include "BioFVM_solvers.h" 
#include "BioFVM_vector.h" 

#include <iostream>
#include <omp.h>

namespace BioFVM{

// do I even need this? 
void diffusion_decay_solver__constant_coefficients_explicit( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay) ... " << std::endl << std::endl;  

		if( M.mesh.uniform_mesh == true )
		{
			std::cout << "Uniform mesh detected! Consider switching to a more efficient method, such as " << std::endl  
			<< "     diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh" << std::endl  
			<< std::endl; 
		}

		precomputations_and_constants_done = true; 
	}

	return; 
}

void diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay, uniform mesh) ... " << std::endl << std::endl;  

		if( M.mesh.uniform_mesh == false )
		{ std::cout << "Error. This code is only supported for uniform meshes." << std::endl; }

		precomputations_and_constants_done = true; 
	}

	return; 
}


void diffusion_decay_solver__constant_coefficients_LOD_3D( Microenvironment& M, double dt )
{
	if( M.mesh.uniform_mesh == false || M.mesh.Cartesian_mesh == false )
	{
		std::cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: other solvers!" << std::endl << std::endl; 
	return; 
	}

	// define constants and pre-computed quantities 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (implicit 3-D LOD with Thomas Algorithm) ... " 
		<< std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );
		
		M.thomas_denomz.resize( M.mesh.z_coordinates.size() , M.zero );
		M.thomas_cz.resize( M.mesh.z_coordinates.size() , M.zero );

		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 
		M.thomas_k_jump = M.thomas_j_jump * M.mesh.y_coordinates.size(); 

		M.thomas_constant1 =  M.diffusion_coefficients; // dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/3)* dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 		
			
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 /= 3.0; // for the LOD splitting of the source 

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 

		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cy.assign( M.mesh.y_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomy.assign( M.mesh.y_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomy[0] = M.thomas_constant3a; 
		M.thomas_denomy[ M.mesh.y_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.y_coordinates.size() == 1 )
		{ M.thomas_denomy[0] = M.one; M.thomas_denomy[0] += M.thomas_constant2; } 

		M.thomas_cy[0] /= M.thomas_denomy[0]; 
		for( int i=1 ; i <= M.mesh.y_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomy[i] , M.thomas_constant1 , M.thomas_cy[i-1] ); 
			M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cz.assign( M.mesh.z_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomz.assign( M.mesh.z_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomz[0] = M.thomas_constant3a; 
		M.thomas_denomz[ M.mesh.z_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.z_coordinates.size() == 1 )
		{ M.thomas_denomz[0] = M.one; M.thomas_denomz[0] += M.thomas_constant2; } 

		M.thomas_cz[0] /= M.thomas_denomz[0]; 
		for( int i=1 ; i <= M.mesh.z_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomz[i] , M.thomas_constant1 , M.thomas_cz[i-1] ); 
			M.thomas_cz[i] /= M.thomas_denomz[i]; // the value at  size-1 is not actually used  
		}	

		M.diffusion_solver_setup_done = true; 
	}

	// x-diffusion 
	
	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( int k=0; k < M.mesh.z_coordinates.size() ; k++ )
	{
		for( int j=0; j < M.mesh.y_coordinates.size() ; j++ )
		{
			// Thomas solver, x-direction

			// remaining part of forward elimination, using pre-computed quantities 
			int n = M.voxel_index(0,j,k);
			(*M.p_density_vectors)[n] /= M.thomas_denomx[0]; 

			for( int i=1; i < M.mesh.x_coordinates.size() ; i++ )
			{
				n = M.voxel_index(i,j,k); 
				axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
				(*M.p_density_vectors)[n] /= M.thomas_denomx[i]; 
			}

			for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
			{
				n = M.voxel_index(i,j,k); 
				naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] ); 
			}

		}
	}

	// y-diffusion 

	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( int k=0; k < M.mesh.z_coordinates.size() ; k++ )
	{
		for( int i=0; i < M.mesh.x_coordinates.size() ; i++ )
		{
   // Thomas solver, y-direction

	// remaining part of forward elimination, using pre-computed quantities 

	int n = M.voxel_index(i,0,k);
	(*M.p_density_vectors)[n] /= M.thomas_denomy[0]; 

	for( int j=1; j < M.mesh.y_coordinates.size() ; j++ )
	{
		n = M.voxel_index(i,j,k); 
		axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_j_jump] ); 
		(*M.p_density_vectors)[n] /= M.thomas_denomy[j]; 
	}

	// back substitution 
	// n = voxel_index( mesh.x_coordinates.size()-2 ,j,k); 

	for( int j = M.mesh.y_coordinates.size()-2 ; j >= 0 ; j-- )
	{
		n = M.voxel_index(i,j,k); 
		naxpy( &(*M.p_density_vectors)[n] , M.thomas_cy[j] , (*M.p_density_vectors)[n+M.thomas_j_jump] ); 
	}

  }
 }

 // z-diffusion 

	M.apply_dirichlet_conditions();
 #pragma omp parallel for 
 for( int j=0; j < M.mesh.y_coordinates.size() ; j++ )
 {
	 
  for( int i=0; i < M.mesh.x_coordinates.size() ; i++ )
  {
   // Thomas solver, y-direction

	// remaining part of forward elimination, using pre-computed quantities 

	int n = M.voxel_index(i,j,0);
	(*M.p_density_vectors)[n] /= M.thomas_denomz[0]; 

	// should be an empty loop if mesh.z_coordinates.size() < 2  
	for( int k=1; k < M.mesh.z_coordinates.size() ; k++ )
	{
		n = M.voxel_index(i,j,k); 
		axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_k_jump] ); 
		(*M.p_density_vectors)[n] /= M.thomas_denomz[k]; 
	}

	// back substitution 

	// should be an empty loop if mesh.z_coordinates.size() < 2 
	for( int k = M.mesh.z_coordinates.size()-2 ; k >= 0 ; k-- )
	{
		n = M.voxel_index(i,j,k); 
		naxpy( &(*M.p_density_vectors)[n] , M.thomas_cz[k] , (*M.p_density_vectors)[n+M.thomas_k_jump] ); 
		// n -= i_jump; 
	}
  }
 }

 return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_2D( Microenvironment& M, double dt )
{
	if( M.mesh.uniform_mesh == false )
	{
		std::cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}

	
	// constants for the linear solver (Thomas algorithm) 

	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );

		
		// define constants and pre-computed quantities 

		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 

		M.thomas_constant1 =  M.diffusion_coefficients; //   dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/2)*dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 
		
		
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 *= 0.5; // for splitting via LOD

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 
		
		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cy.assign( M.mesh.y_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomy.assign( M.mesh.y_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomy[0] = M.thomas_constant3a; 
		M.thomas_denomy[ M.mesh.y_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.y_coordinates.size() == 1 )
		{ M.thomas_denomy[0] = M.one; M.thomas_denomy[0] += M.thomas_constant2; } 


		M.thomas_cy[0] /= M.thomas_denomy[0]; 
		for( int i=1 ; i <= M.mesh.y_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomy[i] , M.thomas_constant1 , M.thomas_cy[i-1] ); 
			M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used  
		}

		M.diffusion_solver_setup_done = true; 
	}

	// set the pointers

	M.apply_dirichlet_conditions();
	// x-diffusion 
	#pragma omp parallel for 
	for( int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		// Thomas solver, x-direction

		// remaining part of forward elimination, using pre-computed quantities 
		int n = M.voxel_index(0,j,0);
		(*M.p_density_vectors)[n] /= M.thomas_denomx[0]; 

		n += M.thomas_i_jump; 
		for( int i=1; i < M.mesh.x_coordinates.size() ; i++ )
		{
			axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
			(*M.p_density_vectors)[n] /= M.thomas_denomx[i]; 
			n += M.thomas_i_jump; 
		}

		// back substitution 
		n = M.voxel_index( M.mesh.x_coordinates.size()-2 ,j,0); 

		for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
		{
			naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] ); 
			n -= M.thomas_i_jump; 
		}
	}

	// y-diffusion 

	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( int i=0; i < M.mesh.x_coordinates.size() ; i++ )
	{
		// Thomas solver, y-direction

		// remaining part of forward elimination, using pre-computed quantities 

		int n = M.voxel_index(i,0,0);
		(*M.p_density_vectors)[n] /= M.thomas_denomy[0]; 

		n += M.thomas_j_jump; 
		for( int j=1; j < M.mesh.y_coordinates.size() ; j++ )
		{
			axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_j_jump] ); 
			(*M.p_density_vectors)[n] /= M.thomas_denomy[j]; 
			n += M.thomas_j_jump; 
		}

		// back substitution 
		n = M.voxel_index( i,M.mesh.y_coordinates.size()-2, 0); 

		for( int j = M.mesh.y_coordinates.size()-2 ; j >= 0 ; j-- )
		{
			naxpy( &(*M.p_density_vectors)[n] , M.thomas_cy[j] , (*M.p_density_vectors)[n+M.thomas_j_jump] ); 
			n -= M.thomas_j_jump; 
		}
	}

	return; 
}

void diffusion_decay_explicit_uniform_rates( Microenvironment& M, double dt )
{
	using std::vector; 
	using std::cout; 
	using std::endl; 

	static int n_jump_i = 1; 
	static int n_jump_j = M.mesh.x_coordinates.size(); 
	static int n_jump_k = M.mesh.x_coordinates.size() * M.mesh.y_coordinates.size(); 

	
	if( !M.diffusion_solver_setup_done )
	{	
		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 
		M.thomas_k_jump = M.thomas_j_jump * M.mesh.y_coordinates.size(); 
	
		M.diffusion_solver_setup_done = true; 
	}
	
	if( M.mesh.uniform_mesh == false )
	{
		cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: something else" << endl << endl; 
		return; 
	}



	// double buffering to reduce memory copy / allocation overhead 

	static vector< vector<double> >* pNew = &(M.temporary_density_vectors1);
	static vector< vector<double> >* pOld = &(M.temporary_density_vectors2);

	// swap the buffers 

	vector< vector<double> >* pTemp = pNew; 
	pNew = pOld; 
	pOld = pTemp; 
	M.p_density_vectors = pNew; 

	static bool reaction_diffusion_shortcuts_are_set = false; 

	static vector<double> constant1 = (1.0 / ( M.mesh.dx * M.mesh.dx )) * M.diffusion_coefficients; 
	static vector<double> constant2 = dt * constant1; 
	static vector<double> constant3 = M.one + dt * M.decay_rates;

	static vector<double> constant4 = M.one - dt * M.decay_rates;

	#pragma omp parallel for
	for( int i=0; i < (*(M.p_density_vectors)).size() ; i++ )
	{
		int number_of_neighbors = M.mesh.connected_voxel_indices[i].size(); 

		double d1 = -1.0 * number_of_neighbors; 

		(*pNew)[i] = (*pOld)[i];  
		(*pNew)[i] *= constant4; 

		for( int j=0; j < number_of_neighbors ; j++ )
		{
			axpy( &(*pNew)[i], constant2, (*pOld)[  M.mesh.connected_voxel_indices[i][j] ] ); 
		}
		vector<double> temp = constant2; 
		temp *= d1; 
		axpy( &(*pNew)[i] , temp , (*pOld)[i] ); 
	}

	return; 
}

};
