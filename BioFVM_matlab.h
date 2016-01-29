/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.0.4) [1]        #
#                                                                           #
# [1] A. Ghaffarizaeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient  #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics, 2015 (accepted).                                       #
#############################################################################
#                                                                           #
#    Bioinformatics, 2015. DOI: 10.1093/bioinformatics/btv730               #
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
#include <cstring>
#include <vector>
#include <iostream>

#include <ctime>
#include <string>
#include <cstring>
#include <cstdio>

#ifndef __BioFVM_matlab_h__
#define __BioFVM_matlab_h__

namespace BioFVM{

// #include "Matrix.h"

// Paul Macklin wrote these based on documentation on the web. 
// So far, matlab v4 is (partially) supported. Exceptions:
//   no complex matrices
//   no sparse matrices 
//   no text matrices

// To save in matlab and make it compatible, make sure you use: 
//   save -v4 <filename> <variable_name>

struct named_vector_data{
std::vector<std::string> names; 
std::vector< std::vector<double> > data; 
};

std::vector< std::vector<double> > read_matlab( std::string filename );
named_vector_data read_matlab_with_names( std::string filename );

bool write_matlab( std::vector< std::vector<double> >& input , std::string filename );
bool write_matlab( std::vector< std::vector<double> >& input , std::string filename , std::vector<std::string>& names );

FILE* write_matlab_header( int rows, int cols, std::string filename, std::string variable_name );  

// output: FILE pointer, and overwrites rows, cols so you know the size 
FILE* read_matlab_header( int* rows, int* cols , std::string filename ); 

};

#endif 
