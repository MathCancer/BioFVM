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

#include "BioFVM_agent_container.h"
#include "BioFVM_basic_agent.h"
#include "BioFVM_vector.h"


namespace BioFVM{


Agent_Container::Agent_Container(){}

void Agent_Container::initialize(int num_voxels){}

void Agent_Container::register_agent( Basic_Agent* agent ){}

void Agent_Container::remove_agent(Basic_Agent* agent ){}
void Agent_Container::add_agent_to_outer_voxel(Basic_Agent* agent){}
void Agent_Container::remove_agent_from_voxel(Basic_Agent* agent, int voxel_index){}
void Agent_Container::add_agent_to_voxel(Basic_Agent* agent, int voxel_index){}
void Agent_Container::update_all_cells(double dt){}

};