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
#include "BioFVM_utilities.h"

namespace BioFVM{
/*
std::string BioFVM_Version; 
std::string BioFVM_URL; 
*/

clock_t program_start_time;
clock_t program_end_time;
clock_t tic_time;
clock_t toc_time;

double clock_t_per_second = (double) CLOCKS_PER_SEC; 
double total_tictoc_time = 0.0; 

void TIC(void)
{ tic_time = clock(); }

void TOC(void)
{
	toc_time = clock(); 
	total_tictoc_time += stopwatch_value(); 
}

void RUNTIME_TIC(void)
{ program_start_time = clock(); }

void RUNTIME_TOC(void)
{ program_end_time = clock(); }

double stopwatch_value(void)
{ return (toc_time - tic_time)/clock_t_per_second; }

double runtime_stopwatch_value(void)
{ return (program_end_time - program_start_time)/clock_t_per_second; }

void display_stopwatch_value( std::ostream& os, double dIn )
{
	int nDays = (int) floor( (double) (dIn / (60.0*60.0*24.0)) );
	int nHours = (int) floor( (double) ( (dIn - nDays*60*60*24) / (60.0*60.0)) );
	int nMinutes = (int) floor( (double) ( (dIn - nDays*60*60*24 - nHours*60*60 ) / (60.0)) );
	double dSeconds = dIn - nDays*60.0*60.0*24.0 - nHours * 60.0*60.0 - nMinutes * 60.0;

	os << nDays << " days, " << nHours << " hours, " 
	  << nMinutes << " minutes, and " << dSeconds << " seconds ";
	return; 
}

std::string format_stopwatch_value( double dIn)
{
	std::string output; 
	output.resize( 1024 ); 
	int nDays = (int) floor( (double) (dIn / (60.0*60.0*24.0)) );
	int nHours = (int) floor( (double) ( (dIn - nDays*60*60*24) / (60.0*60.0)) );
	int nMinutes = (int) floor( (double) ( (dIn - nDays*60*60*24 - nHours*60*60 ) / (60.0)) );
	double dSeconds = dIn - nDays*60.0*60.0*24.0 - nHours * 60.0*60.0 - nMinutes * 60.0;

	sprintf( (char*) output.c_str(),
	"%d days, %d hours, %d minutes, and %2.4f seconds",
	nDays,nHours,nMinutes,dSeconds);

	return output; 
}

double total_stopwatch_time( void )
{ return total_tictoc_time; }

};