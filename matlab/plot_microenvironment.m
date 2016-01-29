% If you want to put titles on the plots, use this syntax: 
% 
% M = read_microenvironment( 'some_file.mat' ); 
% titles{1} = 'cells'; 
% titles{2} = 'blood vessels';
% titles{3} = 'oxygen'; 
% plot_microenvironment( M , titles );  
% 
% To go without the labels, use: 
% 
% M = read_microenvironment( 'some_file.mat' ); 
% plot_microenvironment( M ); 
% 

function plot_microenvironment( M , titles )
plot_titles = [];
if( nargin > 1 )
    plot_titles = titles;  
end

number_of_fields = length( M.data ); 
number_of_plots = number_of_fields; 

if( nargin == 1 )
    for i=1:number_of_fields
        plot_titles{i} = sprintf('field %d',i); 
    end
end


if( mod(number_of_plots,2) == 1 && number_of_fields > 1 )
    number_of_plots = number_of_plots + 1; 
end

width = ceil( sqrt( number_of_plots )) ;
height = number_of_plots / width ;

[x,y,z] = size( M.data{1} ) ;

mid_index = 1+floor( z/2.0 ) ;

for i=1:number_of_fields
    subplot(height,width,i,'align') ; 
    % figure(1) ; 
    h = contourf( M.X , M.Y, M.data{i}(:,:,mid_index)' , 30 ,  'linecolor', 'none'); axis image;  ;   title( plot_titles{i} ); colorbar; 
    xlabel( 'x' ); 
    ylabel( 'y' ); 
    % set( gca ,'CLim', [0 1] ); 
end

return; 