%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use BioFVM in your project, please cite BioFVM and the version     %
% number, such as below:                                                    %
%                                                                           %
% We solved the diffusion equations using BioFVM (Version 1.0.3) [1]        %
%                                                                           %
% [1] A. Ghaffarizaeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient  %
%    parallelized diffusive transport solver for 3-D biological simulations,%
%    Bioinformatics, 2015. DOI: 10.1093/bioinformatics/btv730               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
% Copyright 2015 Paul Macklin and the BioFVM Project                        %
%                                                                           %
% Licensed under the Apache License, Version 2.0 (the "License");           %
% you may not use this file except in compliance with the License.          %
% You may obtain a copy of the License at                                   %
%                                                                           %
%    http://www.apache.org/licenses/LICENSE-2.0                             %
%                                                                           %
% Unless required by applicable law or agreed to in writing, software       %
% distributed under the License is distributed on an "AS IS" BASIS,         %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  %
% See the License for the specific language governing permissions and       %
% limitations under the License.                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
file_name='output_4.000000_0.000100_10.000000.mat';
x=1;y=2;z=3;

needed_plane=[x,y];

crossesction_index=setdiff([1,2,3], needed_plane);
labels={'x','y','z'};
load(file_name);

for i=5:size(multiscale_microenvironment,1)
    m=multiscale_microenvironment;
    
    temp= unique(sort(m(crossesction_index,:)));
    temp_median= temp(floor(length(temp)/2));
    
    m=m(:,m(crossesction_index,:)==temp_median);
    
    tempx= unique(sort(m(needed_plane(1),:)));
    stepx= abs(tempx(1)-tempx(2));
    minx=tempx(1)-stepx/2;
    maxx=tempx(end)+stepx/2;
    tempy= unique(sort(m(needed_plane(2),:)));
    stepy= abs(tempy(1)-tempy(2));
    miny=tempy(1)-stepy/2;
    maxy=tempy(end)+stepy/2;
    
    num_rows= length(tempx);
    num_cols= length(tempx);
    % scaling the x values to range [1:numrows], the y values to range [1:numcols]
    x_scaled= 1+ floor(num_rows*((m(needed_plane(1),:)-minx)/(maxx-minx)));
    y_scaled= 1+ floor(num_cols*((m(needed_plane(2),:)-miny)/(maxy-miny)));
    
    c1=sparse(x_scaled,y_scaled,m(i,:));
    
    full_matrix=full(c1);
    
    figure
    contourf(full_matrix)
    axis image
    shading flat
    colorbar('FontSize', 14)
    title(['substrate' num2str(i-4)], 'FontSize', 14)
    h = colorbar();
    h1 = get( h , 'ylabel' );
    set( h1 , 'string' , 'concentration' )
    set( h1, 'fontsize', 12 )
    set(gca,'XTick',1:num_rows/2-0.5:num_rows)
    set(gca,'YTick',1:num_cols/2-0.5:num_cols)
    set(gca,'YTickLabel',{num2str(miny),num2str((miny+maxy)/2),num2str(maxy)})
    set(gca,'XTickLabel',{num2str(minx),num2str((minx+maxx)/2),num2str(maxx)})
    xlabel([labels{needed_plane(1)},' (\mum)'])
    ylabel([labels{needed_plane(2)},' (\mum)'])
end
