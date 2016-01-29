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

#ifndef __BioFVM_utilities_h__
#define __BioFVM_utilities_h__

#include <iostream>
#include <ctime>
#include <cmath>
#include <string>

namespace BioFVM{

void TIC(void);
void TOC(void);
void RUNTIME_TIC(void);
void RUNTIME_TOC(void);

double stopwatch_value(void);
double runtime_stopwatch_value(void);

// This returns the total number of seconds in successive TIC() / TOC() calls 
// Useful if you use TOC to "pause" your timer (e.g., to exclude file i/o from benchmark times)
double total_stopwatch_time( void );

void display_stopwatch_value( std::ostream& os , double dIn );
std::string format_stopwatch_value( double dIn);
	
};
 
#endif 
