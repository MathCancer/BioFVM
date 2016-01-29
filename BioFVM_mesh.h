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

#ifndef __BioFVM_mesh_h__
#define __BioFVM_mesh_h__

#include <iostream>
#include <vector> 

#include "BioFVM_matlab.h"

namespace BioFVM{

 /*! \brief Voxels are the basic spatial container for densities, which are networked into meshes. 
  * 
  * Voxels are the basic spatial container for a finite volume method. Voxels are connected to 
  * other voxels into a General_Mesh, here most likely a Cartesian_Mesh. Voxel boundaries are Voxel_Faces. 
  * 
  * A Microenvironment Scale will include a network of Voxels (and Voxel_Faces), 
  * with a vector<double> of densities for each Voxel, along with rate constants, etc. 
  * The Scale may also include a vector<double> of flux coefficients for each Voxel_Face. 
 */
 
class Voxel
{

 private:
	friend std::ostream& operator<<(std::ostream& os, const Voxel& v); 
	/*!< outputs the Voxel to an open ostream 
	 * \param os -- the stream 
	 * \param mv -- the voxel you use this friendly friend operator on
	 * Example: Voxel v; 
	 *          cout << v << endl; 
	*/ 

 public:
	Voxel(); 
	int mesh_index; /*!< voxel's index in a General_Mesh */ 

	double volume; /*!< voxel's volume (cubic spatial units) */ 
	std::vector<double> center; /*!< center of volume */
	bool is_Dirichlet;
	void stream_output_with_units( std::ostream& os , std::string units ) const;
};

class Voxel_Face
{
 private:
	friend std::ostream& operator<<(std::ostream& os , const Voxel_Face& vf ); 
	
 public:
	Voxel_Face(); 
	int mesh_index; 
	
	double surface_area; 
	std::vector<double> center; 
	std::vector<double> outward_normal; 
	std::vector<double> inward_normal; 
	
	void stream_output_with_units( std::ostream& os , std::string units ) const;
};

class General_Mesh
{
 private: 
	friend std::ostream& operator<<(std::ostream& os, const General_Mesh& mesh);  
	
	// this stores the indexing of the voxel faces (connect voxel i to voxel j, face stored at k)
	// only for use in a future release
	// std::unordered_map< int,std::unordered_map<int,int> > voxel_face_index_mapping; 
	
 public:
	General_Mesh();  
	
	// [xmin ymin zmin xmax ymax zmax ]
	std::vector<double> bounding_box; 
	
	std::vector<Voxel> voxels; 
	std::vector<Voxel_Face> voxel_faces; 
	// each voxel[k] has a list of connected voxels -- helpful for some numerical methods 
	std::vector< std::vector<int> > connected_voxel_indices; 
	
	int nearest_voxel_index( std::vector<double>& position );   
	bool is_position_valid(double x, double y, double z);
	/* the following help manage the voxel faces */ 

	// returns the index of the voxel face connecting from voxels[i] to voxels[j] 
	int voxel_face_index( int i, int j ); 
	
	// returns the Voxel_Face connecting voxels[i] to voxels[j] 
	Voxel_Face& voxel_face(int i, int j );   
	// returns the normal vector from voxels[i] to voxels[j] 
	std::vector<double>& outward_normal( int i, int j ); 
	
	/*! This creates a Voxel_Face from voxels[i] to voxels[j], and another from voxels[j] to 
	    voxels[i], both with surface area SA. It also auto-updates connected_voxel_indices[i] 
	    and connected_voxel_indices[j]. */ 
	void connect_voxels(int i,int j, double SA);   
	
	void connect_voxels_faces_only(int i,int j, double SA); 
	void connect_voxels_indices_only(int i,int j, double SA); 
	
	/*! This removes all connections between voxels[i] and voxels[j], and deletes the associated 
	    Voxel_Face(s). */
	void disconnect_voxels(int i, int j); 
	void clear_voxel_face_index_mapping( void );  
	
	bool Cartesian_mesh; 
	bool uniform_mesh; 
	bool regular_mesh;
	bool use_voxel_faces; 
	
	std::string units; 
	
	void display_information( std::ostream& os); 
	
	void write_to_matlab( std::string filename ); 
	void read_from_matlab( std::string filename ); 
};

class Cartesian_Mesh : public General_Mesh
{
 private:
 
 public:
	std::vector<double> x_coordinates; 
	std::vector<double> y_coordinates;
	std::vector<double> z_coordinates; 	
	
	int voxel_index( int i, int j, int k ); 
	std::vector<int> cartesian_indices( int n ); 
	
	double dx;
	double dy;
	double dz; 
	
	double dV; 	
	double dS;

	double dS_xy;
	double dS_yz; 
	double dS_xz;
	
	Cartesian_Mesh(); // done 
	
	Cartesian_Mesh( int , int , int );  
	
	void create_voxel_faces( void ); 

	void resize( int,int,int ); 
	void resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , int x_nodes, int y_nodes, int z_nodes ); 
	void resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz ); 
	void resize_uniform( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx ); 
	
	int nearest_voxel_index( std::vector<double>& position );   
	int nearest_voxel_face_index( std::vector<double>& position );  
	std::vector<int> nearest_cartesian_indices( std::vector<double>& position ); 
	Voxel& nearest_voxel( std::vector<double>& position ); 
	
	void display_information( std::ostream& os ); 
	
	void read_from_matlab( std::string filename ); 
};

class Voronoi_Mesh : public General_Mesh
{
 private:
 
 public:
	void display_information( std::ostream& os); 
};

};

#endif