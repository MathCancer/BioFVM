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

#include "BioFVM_vector.h" 

/* some global BioFVM strings */ 

namespace BioFVM{

/* faster operator overloading. multiplication and division are element-wise (Hadamard) */ 

std::vector<double> operator-( const std::vector<double>& v1 , const std::vector<double>& v2 )
{
 std::vector<double> v = v1;
 for( int i=0; i < v1.size() ; i++ )
 { v[i] -= v2[i]; }
 return v; 
}

std::vector<double> operator+( const std::vector<double>& v1 , const std::vector<double>& v2 )
{
 std::vector<double> v = v1;
 for( int i=0; i < v1.size() ; i++ )
 { v[i] += v2[i]; }
 return v; 
}

std::vector<double> operator*( const std::vector<double>& v1 , const std::vector<double>& v2 )
{
 std::vector<double> v = v1;
 for( int i=0; i < v1.size() ; i++ )
 { v[i] *= v2[i]; }
 return v; 
}

std::vector<double> operator/( const std::vector<double>& v1 , const std::vector<double>& v2 )
{
 std::vector<double> v = v1;
 for( int i=0; i < v1.size() ; i++ )
 { v[i] /= v2[i]; }
 return v; 
}

std::vector<double> operator*( double d , const std::vector<double>& v1 )
{
 std::vector<double> v = v1;
 for( int i=0; i < v1.size() ; i++ )
 { v[i] *= d; }
 return v; 
}

std::vector<double> operator+( double d , const std::vector<double>& v1 )
{
 std::vector<double> v = v1;
 for( int i=0; i < v1.size() ; i++ )
 { v[i] += d; }
 return v; 
}

std::vector<double> operator+( const std::vector<double>& v1 , double d )
{
 std::vector<double> v = v1;
 for( int i=0; i < v1.size() ; i++ )
 { v[i] += d; }
 return v; 
}

std::vector<double> operator-( double d , const std::vector<double>& v1 )
{
 std::vector<double> v = v1;
 for( int i=0; i < v1.size() ; i++ )
 { v[i] = d - v1[i]; }
 return v; 
}

std::vector<double> operator-( const std::vector<double>& v1 , double d  )
{
 std::vector<double> v = v1;
 for( int i=0; i < v1.size() ; i++ )
 { v[i] -= d; }
 return v; 
}

void operator+=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( int i=0; i < v1.size() ; i++ )
 { v1[i] += v2[i]; }
 return; 
}

void operator-=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( int i=0; i < v1.size() ; i++ )
 { v1[i] -= v2[i]; }
 return; 
}

void operator/=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( int i=0; i < v1.size() ; i++ )
 { v1[i] /= v2[i]; }
 return;  
} 

void operator*=( std::vector<double>& v1, const double& a )
{
 for( int i=0; i < v1.size() ; i++ )
 { v1[i] *= a; }
 return; 
}

void operator*=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( int i=0; i < v1.size() ; i++ )
 { v1[i] *= v2[i]; }
 return;  
}

void operator/=( std::vector<double>& v1, const double& a )
{
 for( int i=0; i < v1.size() ; i++ )
 { v1[i] /= a; }
 return;  
}

/* other commonly needed operations on vectors */ 

std::ostream& operator<<(std::ostream& os, const std::vector<double>& v )
{
 if( v.size() == 3 )
 {
  os << "x=\"" << v[0] << "\" y=\"" << v[1] << "\" z=\"" << v[2] << "\"" ; 
  return os; 
 }

 for( int i=0; i < v.size(); i++ )
 { os << v[i] << " " ; }
 return os; 
}

// this one returns a new vector that has been normalized
std::vector<double> normalize( std::vector<double>& v )
{
 std::vector<double> output = v ;

 double norm = 0.0; 
 norm = 0.0; 

 for( int i=0; i < v.size(); i++ )
 { norm += ( v[i]*v[i] ); }
 norm = sqrt( norm ); 

 for( int i=0; i < v.size(); i++ )
 { output[i] /= norm ; }
 return output; 
}

// this one normalizes v
void normalize( std::vector<double>* v )
{
 double norm = 1e-32; 

 for( int i=0; i < v->size(); i++ )
 { norm += ( (*v)[i] * (*v)[i] ); }
 norm = sqrt( norm ); 

 for( int i=0; i < v->size(); i++ )
 { (*v)[i] /=  norm ; }
 return; 
}

double norm_squared( const std::vector<double>& v )
{
 double out = 0.0; 
 for( int i=0 ; i < v.size() ; i++ )
 { out += ( v[i] * v[i] ); }
 return out; 
}

double norm( const std::vector<double>& v )
{
 return sqrt( norm_squared( v ) ); 
}

double maxabs( const std::vector<double>& v )
{
 double out = 0.0; 
 for( int i=0; i < v.size() ; i++ )
 {
  if( fabs( v[i] ) > out )
  { out = v[i]; }
 }
 return out; 
}

double max_abs_difference( const std::vector<double>& v1 , const std::vector<double>& v2 )
{
 double out = 0.0; 
 for( int i=0; i < v1.size() ; i++ )
 {
  if( fabs( v1[i] -v2[i] ) > out )
  { out = fabs( v1[i] - v2[i] ); }
 }
 return out; 
}

std::vector<double> exponentiate( const std::vector<double>& exponent )
{
 std::vector<double> out( exponent.size() , 0.0 );
 
 for( int i=0 ; i < out.size() ; i++ )
 { out[i] = exp( exponent[i] ); }

 return out; 
}
 
void randomize( std::vector<double>* v )
{
 double norm = 1e-32; 
 
 static double d1 = 2.0 / (double) RAND_MAX; 

 for( int i=0; i < v->size(); i++ )
 { (*v)[i] =  -1 + d1 * rand(); }
 
 return; 
}

/* axpy and related BLAS-type operations */ 

void axpy( std::vector<double>* y, double& a , std::vector<double>& x )
{
 for( int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] += a * x[i] ; 
 }
 return ; 
}

void axpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x )
{
 for( int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] += a[i] * x[i] ; 
 }
 return; 
}

void naxpy( std::vector<double>* y, double& a , std::vector<double>& x )
{
 for( int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] -= a * x[i] ; 
 }
 return ; 
}

void naxpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x )
{
 for( int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] -= a[i] * x[i] ; 
 }
 return; 
}

// turn a delimited character array (e.g., csv) into a vector of doubles

void csv_to_vector( const char* buffer , std::vector<double>& vect )
{
	vect.resize(0); 
	int i=0;
	while( i < strlen( buffer )  )
	{
		// churn through delimiters, whitespace, etc. to reach the next numeric term
		while( isdigit( buffer[i] ) == false && buffer[i] != '.' && buffer[i] != '-' && buffer[i] != 'e' && buffer[i] != 'E' )
		{ i++; } 
		char* pEnd; 
		if( i < strlen(buffer) ) // add this extra check in case of a final character, e.g., ']'
		{
			vect.push_back( strtod( buffer+i , &pEnd ) ); 
			i = pEnd - buffer; 
		}
	}			
	return; 
}

char* vector_to_csv( const std::vector<double>& vect )
{ 
	static int datum_size = 16;  // format = %.7e, 1 (sign) + 1 (lead) + 1 (decimal) + 7 (figs) + 2 (e, sign) + 3 (exponent) + 1 (delimiter) = 16
	// this is approximately the same at matlab long for single precision. 
	// If you want better precision, use a binary data format like matlab, or (in the future) HDF 

	char* buffer; 
	buffer = new char[ datum_size * vect.size() ];
	
	int position = 0; 
	for( int j=0; j < vect.size()-1 ; j++ )
	{
		position += sprintf( buffer+position , "%.7e," , vect[j] ); 
	}
	sprintf( buffer + position , "%.7e" , vect[ vect.size()-1 ] ); 
	
	return buffer; 
}

void vector_to_csv_safe( const std::vector<double>& vect , char*& buffer )
{ 
	static int datum_size = 16;  // format = %.7e, 1 (sign) + 1 (lead) + 1 (decimal) + 7 (figs) + 2 (e, sign) + 3 (exponent) + 1 (delimiter) = 16
	// this is approximately the same at matlab long for single precision. 
	// If you want better precision, use a binary data format like matlab, or (in the future) HDF 

	if( buffer )
	{ delete [] buffer; } 
	buffer = new char[ datum_size * vect.size() ];
	std::cout << __LINE__ << std::endl; 
	
	int position = 0; 
	for( int j=0; j < vect.size()-1 ; j++ )
	{
		position += sprintf( buffer+position , "%.7e," , vect[j] ); 
	}
	sprintf( buffer + position , "%.7e" , vect[ vect.size()-1 ] ); 
	return; 
}

void vector_to_csv( const std::vector<double>& vect , char*& buffer )
{ 
	// %.7e is approximately the same at matlab longe for single precision. 
	// If you want better precision, use a binary data format like matlab, or (in the future) HDF 

	int position = 0; 
	for( int j=0; j < vect.size()-1 ; j++ )
	{
		position += sprintf( buffer+position , "%.7e," , vect[j] ); 
	}
	sprintf( buffer + position , "%.7e" , vect[ vect.size()-1 ] ); 
	return; 
}


};

