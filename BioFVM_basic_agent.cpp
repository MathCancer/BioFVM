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

#include "BioFVM_basic_agent.h"
#include "BioFVM_agent_container.h"
#include "BioFVM_vector.h" 

namespace BioFVM{

std::vector<Basic_Agent*> all_basic_agents(0); 

Basic_Agent::Basic_Agent()
{
	//give the agent a unique ID  
	static int max_basic_agent_ID = 0; 
	ID = max_basic_agent_ID; // 
	max_basic_agent_ID++; 
	// initialize position and velocity
	position.assign( 3 , 0.0 ); 
	velocity.assign( 3 , 0.0 );
	previous_velocity.assign( 3 , 0.0 ); 
	// link into the microenvironment, if one is defined 
	secretion_rates= new std::vector<double>(0);
	uptake_rates= new std::vector<double>(0);
	saturation_densities= new std::vector<double>(0);
	extern Microenvironment* default_microenvironment;
	
	register_microenvironment( default_microenvironment ); 

	return;	
}

void Basic_Agent::update_position(double dt){ 
//make sure to update current_voxel_index if you are implementing this function
};
bool Basic_Agent::assign_position(std::vector<double> new_position)
{
	return assign_position(new_position[0], new_position[1], new_position[2]);
}

bool Basic_Agent::assign_position(double x, double y, double z)
{
	// std::cout << __FILE__ << " " << __LINE__ << std::endl; 
	if( !get_microenvironment()->mesh.is_position_valid(x,y,z))
	{	
		// std::cout<<"Error: the new position for agent "<< ID << " is invalid: "<<x<<","<<y<<","<<"z"<<std::endl;
		return false;
	}
	position[0]=x;
	position[1]=y;
	position[2]=z;
	current_voxel_index= microenvironment->nearest_voxel_index( position );
	// make sure the agent is not already registered
	get_microenvironment()->agent_container->register_agent(this);
	return true;
}

// Microenvironment& Basic_Agent::current_microenvironment( void )
// { return multiscale_microenvironment->operator[]( selected_microenvironment ); }

// int Basic_Agent::current_microenvironment_index( void )
// { return selected_microenvironment; }

// void Basic_Agent::select_microenvironment( int microenvironment_index )
// { 
// std::cout << __FILE__ << " "<<current_microenvironment_index()<< " "<< __LINE__ << std::endl; 
// selected_microenvironment = microenvironment_index; 
// std::cout << __FILE__ << " "<< __LINE__ << std::endl; 
// } 

// std::vector<double>& Basic_Agent::secretion_rates( int index )
// { return all_secretion_rates[index]; }

// std::vector<double>& Basic_Agent::saturation_densities( int index )  
// { return all_saturation_densities[index]; }

// std::vector<double>& Basic_Agent::uptake_rates( int index ) 
// { return all_uptake_rates[index];  }

// std::vector<double>& Basic_Agent::secretion_rates( void )
// { return all_secretion_rates[selected_microenvironment]; }

// std::vector<double>& Basic_Agent::saturation_densities( void )  
// { return all_saturation_densities[selected_microenvironment]; }

// std::vector<double>& Basic_Agent::uptake_rates( void ) 
// { return all_uptake_rates[selected_microenvironment]; }

void Basic_Agent::set_internal_uptake_constants( double dt )
{
	
		// overall form: dp/dt = S*(T-p) - U*p 
		//   p(n+1) - p(n) = dt*S(n)*T(n) - dt*( S(n) + U(n))*p(n+1)
		//   p(n+1)*temp2 =  p(n) + temp1
		//   p(n+1) = (  p(n) + temp1 )/temp2
		//int nearest_voxel= current_voxel_index;
		double internal_constant_to_discretize_the_delta_approximation = dt * volume / ( (microenvironment->voxels(current_voxel_index)).volume ) ; // needs a fix 
		
// before the fix on September 28, 2015. Also, switched on this day 
// from Delta function sources/sinks to volumetric 
/*		
		// temp1 = dt*S*T 
		cell_source_sink_solver_temp1 = *secretion_rates; 
		cell_source_sink_solver_temp1 *= *saturation_densities; 
		cell_source_sink_solver_temp1 *= dt; 

		// temp2 = 1 + dt*( S + U )
		cell_source_sink_solver_temp2.assign( (*secretion_rates).size() , 1.0 ); 
		axpy( &(cell_source_sink_solver_temp2) , dt , *secretion_rates );
		axpy( &(cell_source_sink_solver_temp2) , dt , *uptake_rates );
*/

		// temp1 = dt*(V_cell/V_voxel)*S*T 
		cell_source_sink_solver_temp1.assign( (*secretion_rates).size() , 0.0 ); 
		cell_source_sink_solver_temp1 += *secretion_rates; 
		cell_source_sink_solver_temp1 *= *saturation_densities; 
		cell_source_sink_solver_temp1 *= internal_constant_to_discretize_the_delta_approximation; 

		// temp2 = 1 + dt*(V_cell/V_voxel)*( S + U )
		cell_source_sink_solver_temp2.assign( (*secretion_rates).size() , 1.0 ); 
		axpy( &(cell_source_sink_solver_temp2) , internal_constant_to_discretize_the_delta_approximation , *secretion_rates );
		axpy( &(cell_source_sink_solver_temp2) , internal_constant_to_discretize_the_delta_approximation , *uptake_rates );	
}

void Basic_Agent::register_microenvironment( Microenvironment* microenvironment_in )
{
	microenvironment = microenvironment_in; 	
	secretion_rates->resize( microenvironment->density_vector(0).size() , 0.0 );
	saturation_densities->resize( microenvironment->density_vector(0).size() , 0.0 );
	uptake_rates->resize( microenvironment->density_vector(0).size() , 0.0 );	

	// some solver temporary variables 
	cell_source_sink_solver_temp1.resize( microenvironment->density_vector(0).size() , 0.0 );
	cell_source_sink_solver_temp2.resize( microenvironment->density_vector(0).size() , 1.0 );
	return; 
}

Microenvironment* Basic_Agent::get_microenvironment( void )
{ return microenvironment; }


Basic_Agent::~Basic_Agent()
{
 return; 
}

Basic_Agent* create_basic_agent( void )
{
	// std::cout << __FILE__ << " create_basic_agent start" << __LINE__ << std::endl; 
	Basic_Agent* pNew; 
	pNew = new Basic_Agent;	 
	all_basic_agents.push_back( pNew ); 
	pNew->index=all_basic_agents.size()-1;
	// std::cout << __FILE__ << " create_basic_agent end" << __LINE__ << std::endl; 
	return pNew; 
}

void delete_basic_agent( int index )
{
	// deregister agent in microenvironment
	all_basic_agents[index]->get_microenvironment()->agent_container->remove_agent(all_basic_agents[index]);
	// de-allocate (delete) the Basic_Agent; 
	
	delete all_basic_agents[index]; 

	// next goal: remove this memory address. 

	// performance goal: don't delete in the middle -- very expensive reallocation
	// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)

	// move last item to index location  
	all_basic_agents[ all_basic_agents.size()-1 ]->index=index;
	all_basic_agents[index] = all_basic_agents[ all_basic_agents.size()-1 ];

	// shrink the vector
	all_basic_agents.pop_back();
	
	return; 
}

void delete_basic_agent( Basic_Agent* pDelete )
{
	// First, figure out the index of this agent. This is not efficient. 

	// int delete_index = 0; 
	// while( all_basic_agents[ delete_index ] != pDelete )
	// { delete_index++; }

	delete_basic_agent(pDelete->index);
	return; 
}

int Basic_Agent::nearest_voxel_index( void )
{ //return microenvironment->nearest_voxel_index( position );  
	return current_voxel_index;
}

std::vector<double>& Basic_Agent::nearest_density_vector( void ) 
{ //return microenvironment->nearest_density_vector( position ); 
	return microenvironment->nearest_density_vector( current_voxel_index ); 
}


void Basic_Agent::simulate_secretion_and_uptake( Microenvironment* pS, double dt )
{
	//int j = pS->nearest_voxel_index( position );
	// std::cout<<"nearest voxel index" <<j<< "; temp1[0]="<< cell_source_sink_solver_temp1<<"; temp2[0]="<< cell_source_sink_solver_temp2<<
		// "; pS(j)="<< (*pS)(j)<<std::endl;
	(*pS)(current_voxel_index) += cell_source_sink_solver_temp1; 
	(*pS)(current_voxel_index) /= cell_source_sink_solver_temp2; 

	return; 
}

// void Basic_Agent::simulate_secretion_and_uptake( double dt , int microenvironment_index )
// {
	// int j = nearest_voxel_index( microenvironment_index ); 
	
	// (*multiscale_microenvironment)[microenvironment_index](j) += cell_source_sink_solver_temp1[microenvironment_index]; 
	// (*multiscale_microenvironment)[microenvironment_index](j) /= cell_source_sink_solver_temp2[microenvironment_index]; 
	
	// return; 
// }

// void Basic_Agent::copy_data( Basic_Agent& copy_me )
// {
	// // copy_rate_data( copy_me );
	// // copy_physical_properties( copy_me );
	
	// copy_phenotype(copy_me);
// }
// void Basic_Agent::copy_phenotype( Basic_Agent& copy_me )
// {
	// phenotype= copy_me.phenotype; // check if there is a need for a copy constructor?
// }

// void Basic_Agent::copy_rate_data( Basic_Agent& copy_me )
// {
	// // resize the microevironment correctly 
	// register_microenvironment( copy_me.get_microenvironment() );
	
	// // now copy the rates!
	
		// secretion_rates = copy_me.secretion_rates; 
		// saturation_densities = copy_me.saturation_densities; 
		// uptake_rates = copy_me.uptake_rates; 
	
	
	// microenvironment = copy_me.microenvironment; 
	// return;
// }

// void Basic_Agent::copy_physical_properties( Basic_Agent& copy_me )
// {
	// type=copy_me.type;
	// radius=copy_me.radius;
	// position = copy_me.position; 
	// velocity = copy_me.velocity; 
	// previous_velocity = copy_me.previous_velocity; 
// }

};