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

#ifndef __BioFVM_vector_h__
#define __BioFVM_vector_h__

#include <iostream>
#include <cstdlib>
#include <vector> 
#include <cmath>
#include <cstring>

namespace BioFVM{

/* faster operator overloading. multiplication and division are element-wise (Hadamard) */ 

std::vector<double> operator-( const std::vector<double>& v1 , const std::vector<double>& v2 );
std::vector<double> operator+( const std::vector<double>& v1 , const std::vector<double>& v2 );
std::vector<double> operator*( const std::vector<double>& v1 , const std::vector<double>& v2 );
std::vector<double> operator/( const std::vector<double>& v1 , const std::vector<double>& v2 );

std::vector<double> operator*( double d , const std::vector<double>& v1 );
std::vector<double> operator+( double d , const std::vector<double>& v1 ); 
std::vector<double> operator+( const std::vector<double>& v1 , double d );
std::vector<double> operator-( double d , const std::vector<double>& v1 );
std::vector<double> operator-( const std::vector<double>& v1 , double d  ); 

void operator+=( std::vector<double>& v1, const std::vector<double>& v2 ); 
void operator-=( std::vector<double>& v1, const std::vector<double>& v2 ); 
void operator/=( std::vector<double>& v1, const std::vector<double>& v2 ); 
void operator*=( std::vector<double>& v1, const double& a );
void operator*=( std::vector<double>& v1, const std::vector<double>& v2 ); 
void operator/=( std::vector<double>& v1, const double& a );

/* other commonly needed operations on vectors */ 

std::ostream& operator<<(std::ostream& os, const std::vector<double>& v ); 

// this one returns a new vector that has been normalized
std::vector<double> normalize( std::vector<double>& v );

// this one normalizes v
void normalize( std::vector<double>* v ); 

double norm_squared( const std::vector<double>& v ); 
double norm( const std::vector<double>& v ); 

double maxabs( const std::vector<double>& v ); 
double max_abs_difference( const std::vector<double>& v1 , const std::vector<double>& v2 ); 

std::vector<double> exponentiate( const std::vector<double>& exponent ); 

// note that the PRNG must be replaced if you are serious about this function
void randomize( std::vector<double>* v ); 

/* axpy and related BLAS-type operations */ 

// y = y + a*x 
void axpy( std::vector<double>* y, double& a , std::vector<double>& x );
// y = y + a.*x
void axpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x ); 

// y = y - a*x 
void naxpy( std::vector<double>* y, double& a , std::vector<double>& x );
// y = y - a.*x
void naxpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x ); 

/* I may cut these from the final version */ 
/* CLEANUP BEFORE RELEASE */ 

//y = y + a.*x  ; y = y ./ d
void axpy_div( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x , std::vector<double>& d );

// y = y + a1.*x1 + a2.*x2  
void double_axpy( std::vector<double>* y, std::vector<double>& a1 , std::vector<double>& a2, std::vector<double>& x1 , std::vector<double>& x2 ); 
// y = y + a.*(x1 + x2)  
void double_axpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x1 , std::vector<double>& x2 ); 
// y = y + a.*(x1 + x2)  ; y = y./d
void double_axpy_div( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x1 , std::vector<double>& x2 , std::vector<double>& d); 
// y = y + a1.*x1 + a2.*x2  ; y = y./d
void double_axpy_div( std::vector<double>* y, std::vector<double>& a1 , std::vector<double>& a2, std::vector<double>& x1 , std::vector<double>& x2 , std::vector<double>& d);

// turn a delimited character array (e.g., csv) into a vector of doubles

void csv_to_vector( const char* buffer , std::vector<double>& vect ); 
char* vector_to_csv( const std::vector<double>& vect );
void vector_to_csv_safe( const std::vector<double>& vect , char*& buffer );
void vector_to_csv( const std::vector<double>& vect , char*& buffer );

};

#endif
