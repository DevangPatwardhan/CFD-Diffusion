
%%% Carteisian structured mesh

% Implement Cell size , grid line/face  
function [nx, ny , gx ,gy , dx , dy ] = meshgen(lx,ly , Nx, Ny ,rx, ry )

%Domain Sizing


% Number of divisions


%Size of the equidistant mesh
dx_size = lx/ Nx;
dy_size = ly/ Ny ;

%Growth/Shrink Factor  using / or * [ Use both rx = ry = 1 for equidistant
%mesh]

%Size of individual lines 
dx = ones(Nx,1) * dx_size ;
dy = ones(Ny,1) * dy_size ;

for i = 1 : Nx 
    dx(i) = dx(i) * ( rx)^i;
end
resize_factorx = sum(dx)/lx;


for i = 1 : Ny
    dy(i) = dy(i) * (ry)^i ;
end
resize_factory = sum(dy)/ly;

%Resizing so that it fits the domain
dx =  dx/resize_factorx;
dy = dy/resize_factory;



%% Grid Lines
gx = zeros(Nx + 1, 1); 
for i = 2:Nx+1
    gx(i) = gx(i-1) + dx(i-1); 
end

gy = zeros(Ny + 1, 1); 

for i = 2:Ny+1
    gy(i) = gy(i-1) + dy(i-1); 
end

%% Node points
nx = zeros(Nx + 2, 1);

for i = 2:Nx+1
    nx(i) = (gx(i-1) + gx(i))/2;  % Taking the midpoint
end
nx(Nx + 2) = gx(Nx + 1 ); % At the end


ny = zeros(Ny + 2, 1);

for i = 2:Ny+1
    ny(i) = (gy(i-1) + gy(i))/2; 
end
ny(Ny + 2) = gy(Ny + 1 );

%node points
[X1, Y1] = meshgrid(nx, ny);


%% Plot the mesh
[x1 , y1 ] = meshgrid( gx , gy ) ;
% figure
% plot(x1, gy, color = 'blue') ; hold on;
% plot(gx, y1' , color = 'blue') ;
% xlabel('Lx')
% ylabel("Ly")
% 
% %Plotting the node points
% plot(X1, ny,'r.');





end



