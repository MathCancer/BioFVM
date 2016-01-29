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

#ifndef __BioFVM_solvers_h__
#define __BioFVM_solvers_h__

#include "BioFVM_microenvironment.h" 

namespace BioFVM{
// /*! diffusion-decay solvers for the equation du/dt = D*Laplacian(u) - lambda*u - U(x)*u + M(X)*(uT-u) */ 

// /*! diffusion-decay solver: 3D LOD implicit (stable method). D and r uniform */  
void diffusion_decay_solver__constant_coefficients_LOD_3D( Microenvironment& M, double dt ); // done
// /*! diffusion-decay solver: 2D LOD implicit (stable method). D and r uniform */  
void diffusion_decay_solver__constant_coefficients_LOD_2D( Microenvironment& M, double dt ); // done

/*! This solves for constant diffusion coefficients on a general mesh using the 
    explicit stepping for the diffusion operator, and implicit stepping for all 
    other terms to increase stability. It is suitable for a general mesh. */ 



// /*! diffusion-decay solver: 3D explicit method -- suitable to a general mesh if necessary */  
void diffusion_decay_explicit_uniform_rates( Microenvironment& M , double dt );  // it exists 

void diffusion_decay_solver__constant_coefficients_explicit( Microenvironment& M, double dt ); 
void diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh( Microenvironment& M, double dt ); 
};

#endif 