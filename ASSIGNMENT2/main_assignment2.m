 clc; clear all ;
% Dimensio
Nx = 80;
Ny = 80;
L = 1;
H = 0.5;
tolerance = 1e-6;
max_time = 10000;


[X_n , Y_n , gx, gy , dx , dy ] = meshgen(L, H ,Nx, Ny , 1/1.1 , 1/1.00);
Y_c = Y_n( 2:end-1); X_c = X_n(2:end-1);
x_coord = X_n ; y_coord = Y_n;
Nx_n = Nx  + 2;
Ny_n = Ny + 2;
%%
area = zeros(length(dx) , length(dy));
% Calculate the area for each grid cell
for i = 1:Nx
    for j = 1:Ny
        area(i,j) = dx(i) * dy(j);
    end
end
totalArea = sum(area(:)); % Sums up correctly


% SOURCE ---------------------------------------------CHECK THISSSSSSS
S_fun= @(x,y,T) -1.5 ;
Sp = 0; 

%K
k = @(x,y) 16 .* (y./H + 1); 


%------------------ Temperature Boundary Conditions ----------------

T1_fun= @(x,y) 15;
%T2_fun = @(x,y) 5 .*( 1 - y ./H) + 15.*sin(pi .*y ./H  );
T2_fun = @(x,y) 5 .*( 1 - y ./H) + 15.*sin(pi .*y ./H  );

T3_fun = @(x,y) 10;
T4_fun = @(x,y) 10; % need not come here


% Initialize temperature field
T = zeros(Nx_n, Ny_n );
T_start = 12.5;                                 %initial guess for Temp
T(:,:) = T_start ;


T(:, 1) = T1_fun(x_coord,0);  % BOTTOM
T(end, :) = T2_fun(L, y_coord); % RIGHT
T(:, end) = T3_fun( x_coord, H  ) ; % TOP
T(1, :) = T4_fun( 0, y_coord ) ; % left boundary



k_ar = zeros(Nx, Ny);

for i=1:Nx
    for j=1:Ny
        k_ar(i, j) = k(X_c(i) ,Y_c(j));
    end
end

time_passed = 1 ;

%______ Gauss Seidel Iteration________________
curr_iter = 1;
tic;
while time_passed < max_time
    R_max = 0;
    T_old = T;
    
    
    for i_n =  2 : Nx_n -1
        for j_n = 2: Ny_n -1
            i = i_n - 1;
            j = j_n - 1;
            
            % Now were ina particular cell (i,j) the center coord of that
            % cell is
            x = X_c(i);
            y = Y_c(j);
            
            % 
            % if i < Nx
            %     kE = k_ar(i +1 , j );
            % else
            %     kE = kP % at the east boundary
            % end
            % kP = k_ar(i,j)
            % 
            % if i> 1 
            %     kW = k_ar(i -1 , j);
            % else
            %     kW = kP; % at the west boundary
            % 
            % end
            % 
            % if j < Ny
            %     kN = k_ar(i,j+1);
            % else
            %     kN = kP; % at the north boundary
            % end
            % 
            % if j > 1
            %     kS = k_ar(i, j -1);
            % else
            %     kS = kP
            % end

            ke = k( gx(i+1) , y  );
             kw = k(gx(i) , y );
              ks = k(x , gy(j) );
               kn = k(x , gy(j+1));
           
            dxe = X_n(i+2) - X_n(i+1);
            dxw = X_n(i+1) - X_n(i);
            dyn = Y_n(j+2) - Y_n(j+1);
            dys = Y_n(j+1) - Y_n(j);
            % 
            % 
            aE = ke * dy(j) / dxe  ;
            aW = kw * dy(j) / dxw;
            aN = kn * dx(i) / dyn;
            aS = ks * dx(i) / dys;

            
            T_cell_old = T_old(i, j);
            cell_area = area(i , j);
            S = S_fun( x,y , T_cell_old);

            aP = aE + aW + aN + aS - Sp;  % _____________________CHECK HERE SP
            
            RHS = aE * T(i_n + 1 ,j_n ) + aW * T(i_n -1 , j_n) + aN * T(i_n , j_n+1 ) + aS* T(i_n,j_n-1) + S* cell_area ;

            T_cell_new = RHS / aP;
            T(i_n, j_n) = T_cell_new; % Update the temperature field
            %%%%%%%%%%%%%%%%%%%%%%%% CHECK HERE>>> OLD OR NEW$$$$$$$$$$$
            cell_residual = aE * T_old(i_n + 1 ,j_n ) + aW * T_old(i_n -1 , j_n) + aN * T_old(i_n , j_n+1 ) + aS* T_old(i_n,j_n-1) + S* cell_area - aP*T_old(i_n,j_n);
            R_max = max(R_max, abs(cell_residual)); % Update the maximum residual
            if mod(curr_iter,100) == 0
            fprintf("Residual %e \n",R_max);
            end


            time_passed = toc;
        end
    end

    %updating the boundary conditoins again to be force it
    T(:, 1) = T1_fun(x_coord,0);  % BOTTOM
    % T(end, :) = T2_fun(L, y_coord); % RIGHT
    
    for j = 1:Ny_n
 
        q_B4 = +5000;                  % prescribed heat flux
        T(end, j) = T(end-1, j) - q_B4*( X_n(end) - X_n(end-1)  )/ke;                     %%%%%%%%%%%%%%%%%%%%%%%YOFKSHIFKJBNSDIKJSBIKVSNFVIKSRHU
    end


    T(:, end) = T3_fun( x_coord, H  ) ; % TOP
    T(1, :)= T(2,:);


Residual(curr_iter) = R_max;
% Check for convergence
    if R_max < tolerance
        fprintf("Convergence achieved after %d iteratioons with residual %e \n",curr_iter , R_max);
     
        break; % Exit the loop if the maximum residual is within the tolerance
    end
    
    if time_passed > max_time
        fprintf("Solver is too slow , Residual reached %e  \n exiting...\n",R_max);
    end



    curr_iter = curr_iter + 1; % Increment the iteration count

end



%%
[X,Y]= meshgrid(X_c,Y_c);

%-------------------------
% POST-PROCESSING & PLOTS
%-------------------------
% Extract interior CV solution (size N_CV_x  x N_CV_y)
T_interior = T(2:Nx_n-1, 2:Ny_n-1); % size (N_CV_x x N_CV_y)
% For plotting we used meshgrid(X_c, Y_c) which is (N_CV_y x N_CV_x)
% So transpose interior for contourf
figure(1);
contourf(X, Y, T_interior', 25);
colorbar;
axis equal tight;
title('Temperature Contour (T)');
xlabel('x');
ylabel('y');

% Convergence plot
figure;
semilogy(1:length(Residual), Residual);
grid on;
xlabel('Iteration');
ylabel('Max Abs Residual');
title(sprintf('Convergence (TOL = %g)', tolerance));

% Compute gradients at CV centers (central difference using T array)
dTdx = zeros(Nx, Ny);
dTdy = zeros(Nx, Ny);

for i_cv = 1:Nx
    for j_cv = 1:Ny
        i_T = i_cv + 1; j_T = j_cv + 1; % index in T
        % central differences safe because we have ghost/boundary nodes
        dTdx(i_cv, j_cv) = (T(i_T+1, j_T) - T(i_T-1, j_T)) / (2 * dx(i_cv));
        dTdy(i_cv, j_cv) = (T(i_T, j_T+1) - T(i_T, j_T-1)) / (2 * dy(j_cv));
    end
end

%  k at  centers
k_C = zeros(Nx, Ny);
for i = 1:Nx
    for j = 1:Ny
        k_C(i, j) = k(X_c(i), Y_c(j));
    end
end

qx = - (k_C .* dTdx); 
qy = - (k_C .* dTdy);

% For quiver use X (N_CV_y x N_CV_x) and qx' (N_CV_y x N_CV_x)
figure(3);
contourf(X, Y, T_interior', 25);
colorbar;
quiver(X, Y, qx', qy'); hold on;

axis equal tight;
title('\bfHeat Flux Vector Plot (\it\bfq\rm\bf)');
xlabel('x');
ylabel('y');

% Print some useful info
T_max = max(T_interior(:));
T_min = min(T_interior(:));
fprintf('T_max = %.6f, T_min = %.6f\n', T_max, T_min);

% Optionally show a few cross-sections
figure(4);
plot(Y_c, T_interior(round(Ny/2), :), '-o'); % centerline vs y
xlabel('x');
ylabel('T (centerline y)');
title('Centerline Temperature vs x');

% End of script


%%

save(sprintf('T_ %d X %d.mat', Nx,Ny),"T","X_n","Y_n","qx",'qy');

% Optionally show a few cross-sections
%figure(5);
plot(X_n, T(:,round(Ny/2) ), '-o');hold on; % centerline vs x
xlabel('x');
ylabel('T (centerline y)');
title('Centerline Temperature vs x');
