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
// #include "BioFVM_strings.h" 

namespace BioFVM{
std::string BioFVM_Version = "1.0.4"; 
std::string BioFVM_URL = "http://BioFVM.MathCancer.org"; 	
	
bool biofvm_dom_initialized = false; 
pugi::xml_document biofvm_doc; 

std::string MultiCellDS_version_string = "0.3.0 beta"; 
std::string MultiCellDS_experimental_snapshot_type_string = "snapshot/experimental";  
int MultiCellDS_experimental_snapshot_code = 3; 

std::string MultiCellDS_simulation_snapshot_type_string = "snapshot/simulation"; 
int MultiCellDS_simulation_snapshot_code = 2; 

std::string MultiCellDS_digital_cell_line_type_string = "digital_cell_line";
int MultiCellDS_digital_cell_line_code = 1; 

std::string MultiCellDS_generic_type_string = "generic";
int MultiCellDS_generic_code = 0; 

bool biofvm_multicellds_default_metadata_initialized = false; 
MultiCellDS_Metadata default_metadata; 

bool save_mesh_as_matlab = false; 
bool save_density_data_as_matlab = true; 

MultiCellDS_Metadata::MultiCellDS_Metadata()
{
	MultiCellDS_type = MultiCellDS_generic_type_string;  
	
	user = "Unspecified"; 
	program = "BioFVM"; 
	program_version = BioFVM_Version; 
	
	bounding_box.resize( 6, 0.0 ); 
	current_time = 0;  
	
	spatial_units = "microns";
	time_units = "min";
	runtime_time_units = "seconds";
	
	current_runtime = 0; 
}


void save_MultiCellDS_snapshot( std::string filename_base , Microenvironment& M ); 

void set_save_biovm_mesh_as_matlab( bool newvalue )
{ save_mesh_as_matlab = newvalue; }

void set_save_biovm_data_as_matlab( bool newvalue )
{ save_density_data_as_matlab = newvalue; }

void read_microenvironment_from_MultiCellDS_xml( Microenvironment& M_destination , std::string filename )
{
	std::cout << "Reading data from file " << filename << " ... " ; 
	pugi::xml_document doc; 
	pugi::xml_parse_result result = doc.load_file( filename.c_str()  );
	
	read_microenvironment_from_MultiCellDS_xml( M_destination , doc ); 
}

void read_microenvironment_from_MultiCellDS_xml( Microenvironment& M_destination , pugi::xml_document& xml_dom )
{
	// find the first microenvironment 
	
	using namespace pugi;
	xml_node root = xml_dom.child("MultiCellDS");
	root = root.child( "microenvironment");
	root = root.child("scale");

	// read all the microenvironments 
	
	int microenvironment_index = -1; 
	while( root )
	{
		M_destination.name= root.attribute("name").value(); 
		
		// read the mesh 
		bool cartesian = true; 
		
		xml_node node = root.child( "mesh" );
		if( strcmp(  node.attribute( "type" ).value()  , "Cartesian" ) != 0 && strcmp(  node.attribute( "type" ).value()  , "cartesian" ) != 0  ) 
		{ cartesian = false; }

		M_destination.mesh.units = node.attribute("units").value(); 
		M_destination.spatial_units = node.attribute("units").value(); 

		// if the dataset doesn't specify uniform or regular at all, assume it's very simple (regular Cartesian)
		if( node.attribute("uniform" ) )
		{ M_destination.mesh.uniform_mesh = node.attribute("uniform").as_bool(); }
		else
		{ M_destination.mesh.uniform_mesh = true; std::cout << __LINE__ << std::endl; } 
		
		if( node.attribute("regular" ) )
		{ M_destination.mesh.regular_mesh = node.attribute("regular").as_bool(); }
		else
		{ M_destination.mesh.regular_mesh = true; }

		// get the bounding box 
		node = node.child( "bounding_box");
		
		int i=0; 
		
		csv_to_vector( node.text().get() , M_destination.mesh.bounding_box ); 
		
		// if Cartesian, try to get the mesh just by reading the x, y, z coordinates 
		if( cartesian == true )
		{
			bool read_coordinates = true; 
			
			// read the x coordinates 
			node = node.parent();
			node = node.child("x_coordinates"); 
			M_destination.mesh.x_coordinates.clear();
			i=0;
			
			csv_to_vector( node.text().get() , M_destination.mesh.x_coordinates ); 
			
			if( M_destination.mesh.x_coordinates.size() > 1 )
			{ M_destination.mesh.dx = M_destination.mesh.x_coordinates[1] - M_destination.mesh.x_coordinates[0]; }
			else
			{ M_destination.mesh.dx = 1.0; }
		
			// read the y coordinates 
			node = node.parent();
			node = node.child("y_coordinates"); 
			M_destination.mesh.y_coordinates.clear();
			i=0;
			
			csv_to_vector( node.text().get() , M_destination.mesh.y_coordinates );
			
			if( M_destination.mesh.y_coordinates.size() > 1 )
			{ M_destination.mesh.dy = M_destination.mesh.y_coordinates[1] - M_destination.mesh.y_coordinates[0]; }
			else
			{ M_destination.mesh.dy = 1.0; }

			// read the z coordinates 
			node = node.parent();
			node = node.child("z_coordinates"); 
			M_destination.mesh.z_coordinates.clear();
			i=0;
			csv_to_vector( node.text().get() , M_destination.mesh.z_coordinates );
			if( M_destination.mesh.z_coordinates.size() > 1 )
			{ M_destination.mesh.dz = M_destination.mesh.z_coordinates[1] - M_destination.mesh.z_coordinates[0]; }
			else
			{ M_destination.mesh.dz = 1.0; }
			
			M_destination.mesh.dV = M_destination.mesh.dx * M_destination.mesh.dy * M_destination.mesh.dz; 
			M_destination.mesh.dS = M_destination.mesh.dx * M_destination.mesh.dy; 
			M_destination.mesh.dS_xy = M_destination.mesh.dx * M_destination.mesh.dy;
			M_destination.mesh.dS_yz = M_destination.mesh.dy * M_destination.mesh.dz;
			M_destination.mesh.dS_xz = M_destination.mesh.dx * M_destination.mesh.dz;
			
			// now, use this mesh information to properly initialize M_destination.mesh
			
			
			if( M_destination.mesh.x_coordinates.size() < 2 )
			{ M_destination.mesh.dx = M_destination.mesh.bounding_box[3] - M_destination.mesh.bounding_box[0]; } 
			if( M_destination.mesh.y_coordinates.size() < 2 )
			{ M_destination.mesh.dy = M_destination.mesh.bounding_box[4] - M_destination.mesh.bounding_box[1]; } 
			if( M_destination.mesh.z_coordinates.size() < 2 )
			{ M_destination.mesh.dz = M_destination.mesh.bounding_box[5] - M_destination.mesh.bounding_box[2]; } 

			if( M_destination.mesh.regular_mesh || M_destination.mesh.uniform_mesh )
			{
				M_destination.resize_space( M_destination.mesh.bounding_box[0], M_destination.mesh.bounding_box[3], 
					M_destination.mesh.bounding_box[1], M_destination.mesh.bounding_box[4], 
					M_destination.mesh.bounding_box[2], M_destination.mesh.bounding_box[5],
					M_destination.mesh.dx, M_destination.mesh.dy, M_destination.mesh.dz ); 			
			}
			else
			{
				std::cout << "Warning! BioFVM / MultiCellDS cannot currently fully initialize to non-regular Cartesian data!" << std::endl; 
			}
		}
		
		if( cartesian == false || M_destination.mesh.regular_mesh == false )
		{
			// Read in the voxels here and create them.
			// If non-regular Cartesian, create and populate them here, one by one. 
			
			node = node.parent(); // now we're at the mesh level 
			pugi::xml_node node1 = node; // set a "bookmark" at the mesh level 

			node = node.child( "voxels" );  // now we're at the voxels level 
			
			// are the voxels written in matlab format or as xml? 
			pugi::xml_attribute attrib = node.attribute( "type" );
			
			if( strcmp(  node.attribute( "type" ).value()  , "matlab" ) == 0 )
			{
				std::cout << "matlab" << std::endl; 				
				M_destination.mesh.Cartesian_mesh = false; 
				M_destination.mesh.uniform_mesh = false; 
				M_destination.mesh.regular_mesh = false; 
				M_destination.mesh.use_voxel_faces = false; 				
				
				// determine the number of voxels 
				int rows; 
				int columns; 
				FILE* fp = read_matlab_header( &rows, &columns, node.text().get() ); 
				int voxel_count = columns; 
				
				// resize the appropriate data structure 
				M_destination.resize_voxels( voxel_count );	
				

				// read the data directly into the voxels  
				for( int j=0; j < columns ; j++ )
				{
					double temp; 
					// read x, y, z, dV
					fread( (char*) & (M_destination.mesh.voxels[j].center[0])   , sizeof(double) , 1 , fp );
					fread( (char*) & (M_destination.mesh.voxels[j].center[1])   , sizeof(double) , 1 , fp );
					fread( (char*) & (M_destination.mesh.voxels[j].center[2])   , sizeof(double) , 1 , fp );
					fread( (char*) & (M_destination.mesh.voxels[j].volume)   , sizeof(double) , 1 , fp );
				} 
				fclose( fp );				
			}
			else
			{
				// Need to set the mesh to non-Cartesian. 
				// We're in very, very basic mode here. 
				
				M_destination.mesh.Cartesian_mesh = false; 
				M_destination.mesh.uniform_mesh = false; 
				M_destination.mesh.regular_mesh = false; 
				M_destination.mesh.use_voxel_faces = false; 
				
				// first, figure out how many voxels. 
				node = node.child( "voxel" ); 
				std::cout << node.name() << std::endl; 
				int voxel_count = 0; 
				while( node )
				{
					voxel_count++; 
					node = node.next_sibling( "voxel" ); 
				}
				
				// now, resize the data structures 
				M_destination.resize_voxels( voxel_count ); 
				
				// now, go back and read in the data 
				node = root; // microenvironment; 
				node = node.child( "mesh" ); 
				node = node.child( "voxels" );
				node = node.child( "voxel"); 
				
				int voxel_index = 0; 
				while( node ) 
				{
					M_destination.mesh.voxels[voxel_index].mesh_index = node.attribute( "ID" ).as_int(); 
					
					// now, get the coordinates 
					node = node.child("center"); 
					csv_to_vector( node.first_child().value() , M_destination.mesh.voxels[voxel_index].center ); 
					node = node.parent(); 
					
					// now, get the volume
					node = node.child( "volume"); 
					M_destination.mesh.voxels[voxel_index].volume = strtod( node.first_child().value() , NULL );  
					node = node.parent(); 

					voxel_index++; 
					node = node.next_sibling( "voxel" ); 
				}	
			}
			
			node = node1; // now we're at the mesh level 
		}
		else
		{
			node = node.parent(); // now we're at the mesh level 
		}	
		
		// if not Cartesian, or if we couldn't read the coordinates, read all the voxels manually 
		
		// after reading the mesh, get the data
		node=node.parent(); // now we're at the microenvironment level 
		node=node.child("densities"); 
		
		// first, read in all the densities names and properties, and resize the data appropriately. 

		node=node.child("density");
		bool added_first_density = false; 
		int substrate_index = 0; 
		while( node )
		{
			if( added_first_density == false )
			{
				M_destination.set_density( 0 , node.attribute("name").value() , node.attribute("units").value() );
				added_first_density = true;
				substrate_index = 0; 
			}
			else
			{
				M_destination.add_density( node.attribute("name").value() , node.attribute("units").value() ); 
				substrate_index++; 
			}
			
			// get the diffusivity 
			node = node.child( "diffusion_constant" ); 
			M_destination.diffusion_coefficients[substrate_index] = strtod(  node.first_child().value() , NULL ); 
			node = node.parent();

			// get the decay rate 
			node = node.child( "decay_rate" ); 
			M_destination.decay_rates[substrate_index] = strtod(  node.first_child().value() , NULL ); 
			node = node.parent(); 
			
			node = node.next_sibling("density");
		}
		
		// lastly, read in all the density data 
		node = root.child( "data" ); 
		
		// read in if stored as matlab 		
		if( strcmp(  node.attribute( "type" ).value()  , "matlab" ) == 0 ) 
		{  
			int rows; 
			int columns; 
			FILE* fp = read_matlab_header( &rows, &columns, node.text().get() ); 			
			int start_row = 0; 
			if( rows > M_destination.number_of_densities() )
			{ start_row = 4; }
			

			// read the data directly into the microenvironment 
			for( int j=0; j < columns ; j++ )
			{
				double temp; 
				// read x,y,z,dV to a temp (throwaway) variable (if start_row == 4)
				for( int i=0; i < start_row ; i++ )
				{ fread( (char*) &temp , sizeof(double) , 1 , fp ); }

				// now, read the actual data 
				for( int i=start_row; i < rows ; i++ )
				{ fread( (char*) &( M_destination.density_vector(j)[i-start_row] ) , sizeof(double) , 1 , fp ); }
			} 
			
			fclose( fp );
		}	
		else
		{
			// attempt to read it in as XML data, voxel by voxel 
			node = node.child( "density_vector" ); 
			for( int j=0 ; j < M_destination.mesh.voxels.size() ; j++ )
			{
				csv_to_vector( node.first_child().value() , M_destination.density_vector(j)  ); 
				if( node.next_sibling( "density_vector" ) ) 
				{ node = node.next_sibling( "density_vector" ); }		
			}
			node = node.parent(); 
		}
		root = root.next_sibling(); 
	}		

	std::cout << "done!" << std::endl; 
	return; 
} 

void save_microenvironment_to_MultiCellDS_xml_pugi( std::string filename_base , Microenvironment& M ) 
{
	std::cout << "Writing data to file " << filename_base << ".xml ... ";

	// if biofvm_doc is already initialized with all the mesh information, just update the metadata, 
	// update the data, and save 
	
	// if the DOM has already been created, then merely update it 
	pugi::xml_node root = biofvm_doc.child( "MultiCellDS");
	int default_microenvironment_index=0;
	if( root )
	{
		// simulation time 
		pugi::xml_node node = root.child( "metadata");
		node = node.child( "current_time" ); 
		char buffer [1024];
		sprintf( buffer , "%f" , default_metadata.current_time ); 
		node.first_child().set_value( buffer ); // default_metadata.current_time ); 

		// current wall time 
		node = node.parent(); 
		node = node.child( "current_runtime" ); 
		sprintf( buffer , "%f" , default_metadata.current_runtime ); 
		node.first_child( ).set_value( buffer ); 
		
		node = root.child( "microenvironment" ); 
		node = node.child( "scale" ); 
		
		
		// and now update the density data 
		
		// iterate over the microenvironments 
		
		node = node.child( "data" ); 
		pugi::xml_node node_microenvironment = node; 
		
		if( save_density_data_as_matlab == true )
		{
			// say where the data are stored, and store them;
			char filename [1024]; 
			sprintf( filename , "%s_microenvironment%d.mat" , filename_base.c_str() , default_microenvironment_index ); 
			M.write_to_matlab( filename ); 
			
			pugi::xml_attribute attrib = node.attribute( "type" ); 
			attrib.set_value( "matlab" ); 
			
			node.first_child().set_value( filename );
			
		}
		else
		{
			pugi::xml_attribute attrib = node.attribute( "type" );  
			attrib.set_value( "xml" ); 
			
			node = node.child( "density_vector" ); 
			
			int datum_size = 16; // enough for sprintf default 6 decimal places + period + 5 leading figs (safety) + delimiter + 2 chars safety = 15
			int data_size = datum_size * M.number_of_densities(); 
			
			char* buffer; 
			buffer = new char [data_size]; 
			for( int j=0 ; j < M.mesh.voxels.size() ; j++ )
			{
				vector_to_csv( M.density_vector(j) , buffer ); 
				
				node.first_child().set_value( buffer ); 
				if( node.next_sibling( "density_vector" ) ) 
				{ node = node.next_sibling( "density_vector" ); }					
			}
			delete [] buffer; 
			node = node.parent(); // back up to level of data 
		}		
		
		node = node.parent(); // back to the level of microenvironment 
		node = node.next_sibling( "scale" ); 
			
	}
	else
	{		
		// create the header 
		
		root = biofvm_doc.append_child( "MultiCellDS" ); 
		pugi::xml_attribute attrib = root.append_attribute( "version" ); 
		attrib.set_value( MultiCellDS_version_string.c_str() ); 
		attrib = root.append_attribute( "type" ); 
		attrib.set_value( MultiCellDS_simulation_snapshot_type_string.c_str() ); 

		// add metadata 
		pugi::xml_node node = root.append_child( "metadata" ); 
		// user 
		node = node.append_child( "user" ); 
		node.append_child( pugi::node_pcdata ).set_value( default_metadata.user.c_str() ); 
		// program
		node = node.parent();
		node = node.append_child( "program" ); 
		attrib = node.append_attribute( "version" ); 
		attrib.set_value( default_metadata.program_version.c_str() ); 
		node.append_child( pugi::node_pcdata ).set_value( default_metadata.program.c_str() ); 
		// bounding_box (may cut in future editions -- already specified in the microenvironments)
		node = node.parent(); 
		node = node.append_child("bounding_box"); 
		attrib = node.append_attribute( "type" ); 
		attrib.set_value( "axis-aligned" ); 
		attrib = node.append_attribute( "units" ); 
		attrib.set_value( default_metadata.spatial_units.c_str() ); 
		
		char* buffer;
		buffer = new char [1024];
		sprintf( buffer , "[%f %f %f %f %f %f]" , default_metadata.bounding_box[0] , default_metadata.bounding_box[1] , default_metadata.bounding_box[2] , 
			default_metadata.bounding_box[3] , default_metadata.bounding_box[4] , default_metadata.bounding_box[5] ); 
		node.append_child( pugi::node_pcdata ).set_value( buffer ); 
		// delete buffer; 

		// simulation time 
		node = node.parent(); 
		node = node.append_child( "current_time" ); 
		attrib = node.append_attribute( "units" ); 
		attrib.set_value( default_metadata.time_units.c_str() ); 
		sprintf( buffer , "%f" , default_metadata.current_time ); 
		node.append_child( pugi::node_pcdata ).set_value( buffer ); 

		// current wall time 
		node = node.parent(); 
		node = node.append_child( "current_runtime" ); 
		attrib = node.append_attribute( "units" ); 
		attrib.set_value( default_metadata.runtime_time_units.c_str() ); 
		sprintf( buffer , "%f" , default_metadata.current_runtime ); 
		node.append_child( pugi::node_pcdata ).set_value( buffer ); 
		
		root = root.append_child( "microenvironment" ); 
		
		// now, add microenvironments 		
			node = root.append_child( "scale" ); 
			attrib = node.append_attribute( "name" );
			attrib.set_value( M.name.c_str() ); 
			
			// add mesh
			
			node = node.append_child( "mesh" ); 

			// information about the mesh 
			attrib = node.append_attribute( "type" ); 
			if( M.mesh.Cartesian_mesh == true )
			{ attrib.set_value( "Cartesian" ); }
			else
			{ attrib.set_value( "general"); }
			attrib = node.append_attribute( "uniform" ); 
			attrib.set_value( M.mesh.uniform_mesh ); 
			attrib = node.append_attribute( "regular" ); 
			attrib.set_value( M.mesh.regular_mesh ); 		
			attrib = node.append_attribute( "units" ); 
			attrib.set_value( M.mesh.units.c_str() );

			// add the bounding box 

			node = node.append_child( "bounding_box" ); 
			attrib = node.append_attribute( "type" ); 
			attrib.set_value( "axis-aligned" ); 
			attrib = node.append_attribute( "units" ); 
			attrib.set_value( M.mesh.units.c_str() ); 
			
			sprintf( buffer , "[%f %f %f %f %f %f]" , M.mesh.bounding_box[0] , M.mesh.bounding_box[1] , M.mesh.bounding_box[2] , 
			M.mesh.bounding_box[3] , M.mesh.bounding_box[4] , M.mesh.bounding_box[5] ); 
			node.append_child( pugi::node_pcdata ).set_value( buffer ); 
			node = node.parent(); 

			// if Cartesian, add the x, y, and z coordinates 
			if( M.mesh.Cartesian_mesh )
			{
				char temp [10240];
				int position = 0; 
				for( int k=0 ; k < M.mesh.x_coordinates.size()-1 ; k++ )
				{ position += sprintf( temp+position, "%f," , M.mesh.x_coordinates[k] ); }
				sprintf( temp+position , "%f" , M.mesh.x_coordinates[ M.mesh.x_coordinates.size()-1] ); 
				node = node.append_child( "x_coordinates" ); 
				node.append_child( pugi::node_pcdata ).set_value( temp ); 
				attrib = node.append_attribute("delimiter");
				attrib.set_value( "," ); 
				
				node = node.parent();
				position = 0; 
				for( int k=0 ; k < M.mesh.y_coordinates.size()-1 ; k++ )
				{ position += sprintf( temp+position, "%f," , M.mesh.y_coordinates[k] ); }
				sprintf( temp+position , "%f" , M.mesh.y_coordinates[ M.mesh.y_coordinates.size()-1] ); 
				node = node.append_child( "y_coordinates" ); 
				node.append_child( pugi::node_pcdata ).set_value( temp ); 
				attrib = node.append_attribute("delimiter");
				attrib.set_value( "," ); 
				
				node = node.parent();
				position = 0; 
				for( int k=0 ; k < M.mesh.z_coordinates.size()-1 ; k++ )
				{ position += sprintf( temp+position, "%f," , M.mesh.z_coordinates[k] ); }
				sprintf( temp+position , "%f" , M.mesh.z_coordinates[ M.mesh.z_coordinates.size()-1] ); 
				node = node.append_child( "z_coordinates" ); 
				node.append_child( pugi::node_pcdata ).set_value( temp ); 
				attrib = node.append_attribute("delimiter");
				attrib.set_value( "," ); 
				node = node.parent(); 
			}
			
			// write out the voxels -- minimal data, even if redundant for cartesian 
			if( save_mesh_as_matlab == false )
			{
				node = node.append_child( "voxels" ); 
				attrib = node.append_attribute("type");
				attrib.set_value( "xml"); 
				char temp [1024]; 
				for( int k=0; k < M.mesh.voxels.size() ; k++ )
				{
					node = node.append_child( "voxel" );
					
					attrib = node.append_attribute( "ID" ); 
					attrib.set_value( M.mesh.voxels[k].mesh_index ); 
					attrib = node.append_attribute( "type" ); 
					attrib.set_value( "cube" ); // allowed: cube or unknown 

					node = node.append_child( "center" );
					attrib = node.append_attribute( "delimiter" );
					attrib.set_value( "," );
					sprintf( temp , "%f,%f,%f" , M.mesh.voxels[k].center[0] , M.mesh.voxels[k].center[1], M.mesh.voxels[k].center[2] );
					node.append_child( pugi::node_pcdata ).set_value( temp ); 
					node = node.parent(); 
					
					node = node.append_child( "volume" );
					sprintf( temp , "%f" , M.mesh.voxels[k].volume );
					node.append_child( pugi::node_pcdata ).set_value( temp ); 
					node = node.parent(); 

					node = node.parent(); 
				}
				node = node.parent(); 
			}
			else
			{
				node = node.append_child( "voxels" ); 
				attrib = node.append_attribute("type");
				attrib.set_value( "matlab"); 
		
				char filename [1024]; 
				sprintf( filename , "%s_mesh%d.mat" , filename_base.c_str() , default_microenvironment_index ); 
				M.mesh.write_to_matlab( filename ); 
				node.append_child( pugi::node_pcdata ).set_value( filename ); 
				
				node = node.parent(); 
			}
			
			
			// now, write out the information on the densities  
			// mental model: list of materials and material properties 

			node = node.parent(); // back to level of "microenvironment"  
			node = node.append_child( "densities" ); 
			
			char temp [1024]; 
			for( int j=0 ; j < M.number_of_densities() ; j++ )
			{
				node = node.append_child( "density" ); 
				attrib = node.append_attribute( "name" ); 
				attrib.set_value( M.density_names[j].c_str() ); 
				attrib = node.append_attribute( "units" ); 
				attrib.set_value( M.density_units[j].c_str() ); 
				
				node = node.append_child( "diffusion_constant" ); 
				sprintf( temp , "%f" , M.diffusion_coefficients[j] ); 
				node.append_child( pugi::node_pcdata ).set_value( temp );
				node = node.parent(); 

				node = node.append_child( "decay_rate" ); 
				sprintf( temp , "%f" , M.decay_rates[j] ); 
				node.append_child( pugi::node_pcdata ).set_value( temp );
				node = node.parent(); 

			   node = node.parent(); 
			}
			
			// now write the actual density data 
			
			node = node.parent(); // back to level of "microenvironment"  
			node = node.append_child( "data" ); 
			
			if( save_density_data_as_matlab == true )
			{
				// say where the data are stored, and store them;
				char filename [1024]; 
				sprintf( filename , "%s_microenvironment%d.mat" , filename_base.c_str() , default_microenvironment_index ); 
				M.write_to_matlab( filename ); 
				
				attrib = node.append_attribute( "type" ); 
				attrib.set_value( "matlab" ); 
				
				node.append_child( pugi::node_pcdata ).set_value( filename );				
			}
			else
			{				
				attrib = node.append_attribute( "type" ); 
				attrib.set_value( "xml" ); 
				int datum_size = 16; // enough for sprintf default 6 decimal places + period + 5 leading figs (safety) + delimiter + 2 chars safety = 15
				int data_size = datum_size * M.number_of_densities(); 
				
				char* buffer; 
				buffer = new char [data_size]; 
				for( int j=0 ; j < M.mesh.voxels.size() ; j++ )
				{
					vector_to_csv( M.density_vector(j) , buffer ); 
					node = node.append_child( "density_vector"); 
					attrib = node.append_attribute( "voxel_ID" ); 
					attrib.set_value( M.mesh.voxels[j].mesh_index ); 
					attrib = node.append_attribute( "delimiter" ); 
					attrib.set_value( "," ); 
					
					node.append_child( pugi::node_pcdata ).set_value( buffer ); 
					node = node.parent(); 
					
				}
				delete [] buffer; 
			}		
	}	
	
	char filename[1024]; 
	sprintf( filename , "%s.xml" , filename_base.c_str() ); 
	biofvm_doc.save_file( filename );
	
	std::cout << "done!" << std::endl; 
	return; 
} 

};