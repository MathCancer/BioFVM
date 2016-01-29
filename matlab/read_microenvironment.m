function out = read_microenvironment( filename )

str = sprintf('Reading microenvironment scale stored in %s.\n\tTo access data field (n) at (i,j,k), use out.data{n}(i,j,k). \n' , filename);
disp(str);  

A = struct2array( load( filename ) );

% figure out X, Y, Z, and size

% out = struct([]); 
 
xmin = A(1,1); 
ymin = A(2,1); 
zmin = A(3,1); 

n = size(A,2);
xmax = A(1,n); 
ymax = A(2,n); 
zmax = A(3,n); 

% figure out number of x nodes  
xnodes = 1; 
while( A(1,xnodes) < xmax - eps )
xnodes = xnodes+1; 
end

out.X = A(1,1:xnodes); 

% figure out number of y nodes 
ynodes = 1; 

while( A(2,ynodes*xnodes) < ymax - eps )
ynodes = ynodes + 1;
end

out.Y = A(2,1:xnodes:xnodes*ynodes);

% figure out number of z nodes 

znodes = 1; 
while( A(3,ynodes*xnodes*znodes) < zmax - eps )
znodes = znodes + 1; 
end

out.Z = A(3,1:xnodes*ynodes:xnodes*ynodes*znodes);

% read in data 

temp = zeros(xnodes,ynodes,znodes);
% out.data = [] ; 
% out.data = temp; 
% out.data(1)  = {temp}; 

% out.data = struct([]);


% to access x(i),y(k),z(k) of data field n: 
% out.data{n}(i,j,k) 

for n=5:size(A,1)
    
    
    % fill in the data 
    m = 1; 
    for k=1:znodes
        for j=1:ynodes
            for i=1:xnodes 
                temp(i,j,k) = A(n,m); 
                m = m+1; 
            end
        end
    end

    % allocate the spot 
	out.data(n-4) = {temp}; % {temp}; 
    
    
end









return ;