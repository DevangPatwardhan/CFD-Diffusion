 clear all ;
% Dimension
Nx = 40;
Ny = 40;
rx = 1/1.05;
ry = 1.0;
L = 1;
H = 0.5;
tolerance = 1e-6;
max_time = 10000;


[X_n , Y_n , gx, gy , dx , dy ] = meshgen(L, H ,Nx, Ny , rx , ry);
Y_c = Y_n( 2:end-1); X_c = X_n(2:end-1);
x_coord = X_n ; y_coord = Y_n;
Nx_n = Nx  + 2;
Ny_n = Ny + 2;


[X1, Y1] = meshgrid(X_n, Y_n);
[x1 , y1 ] = meshgrid( gx , gy ) ;
figure(8);
plot(x1, gy, color = 'blue') ; hold on;
plot(gx, y1' , color = 'blue') ;
xlabel('Lx')
ylabel("Ly")
%Plotting the node points
plot(X1, Y_n,'r.');



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
T2_fun = @(x,y) -15 ;

T3_fun = @(x,y) 10;
T4_fun = @(x,y) 10; % need not come here


% Initialize temperature field
T = zeros(Nx_n, Ny_n );
T_start = 12.5;                                 %initial guess for Temp
T(:,:) = T_start ;


T(:, 1) = T1_fun(x_coord,0);  % BOTTOM
%T(end, :) = T2_fun(L, y_coord); % RIGHT
%T(:, end) = T3_fun( x_coord, H  ) ; % TOP
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

F = 1.5;
tic;
while time_passed(end) < max_time
    R_max = 0;
    T_old = T;
    R_sum = 0;
    
    
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
            R_sum = R_sum + abs(cell_residual);
            R_max = max(R_max, abs(cell_residual)); 


        end
    end
    
    
    time_passed = toc;
    %updating the boundary conditoins again to be force it
    T(:, 1) = T1_fun(x_coord,0);  % BOTTOM
    % T(end, :) = T2_fun(L, y_coord); % RIGHT
    
    for j = 1:Ny_n

        q_B2 = +5000;                  % prescribed heat flux
        T(end, j) = T(end-1, j) - q_B2*( X_n(end) - X_n(end-1)  )/k(gx(end), Y_n(j) );                     %%%%%%%%%%%%%%%%%%%%%%%YOFKSHIFKJBNSDIKJSBIKVSNFVIKSRHU
    end

    % T(end , :) = T2_fun(L,y_coord   );



    T(:, end) = T3_fun( x_coord, H  ) ; % TOP

    for i = 1:Nx_n

        q_B3 = +5000;                  % prescribed heat flux
        T(i, end) = T(i,end-1) - q_B2*( Y_n(end) - Y_n(end-1)  )/k(X_n(i), gy(end) );                     %%%%%%%%%%%%%%%%%%%%%%%YOFKSHIFKJBNSDIKJSBIKVSNFVIKSRHU
    end

    T(1, :)= T(2,:);  % Left

    if mod(curr_iter,1000) == 0
        %fprintf("Residual %e  , Iteration number : %d , Time Elapsed : %f \n", R_sum,curr_iter, time_passed);
    end
Residual(curr_iter) = R_sum;
% Check for convergence
    if R_sum/F < tolerance
        fprintf("For %d x %d mesh , grx-%.2f Convergence achieved after %d iteratioons with residual %e . Time Taken :%f  \n", Nx,Ny,rx,curr_iter , R_sum/F ,time_passed);
     
        break; 
    end
    
    if time_passed(end) > max_time
        fprintf("Solver is  slow , Residual reached %e  \n exiting...\n with time %f ",R_max , time_passed);
    end

    curr_iter = curr_iter + 1; % Increment the iteration count

end



%%
[X,Y]= meshgrid(X_n,Y_n);

T_interior = T(2:Nx_n-1, 2:Ny_n-1); % size (N_CV_x x N_CV_y)

figure(1);
contourf(X, Y, T', 25);
colorbar;
axis equal tight;
title('Temperature Contour   (T)');
xlabel('x');
ylabel('y');

% Convergence plot
figure(2);
semilogy(1:length(Residual), Residual);
grid on;
xlabel('Iteration  ');
ylabel('Residual ');
title(sprintf('Convergence   (tolerance = %g)', tolerance));

% Compute gradients at CV centers (central difference using T array)
dTdx = zeros(Nx_n, Ny_n);
dTdy = zeros(Nx_n, Ny_n);

for i_cv = 1:Nx
    for j_cv = 1:Ny
        i_T = i_cv + 1; j_T = j_cv + 1; 
        dTdx(i_cv, j_cv) = (T(i_T+1, j_T) - T(i_T-1, j_T)) / (2 * dx(i_cv));
        dTdy(i_cv, j_cv) = (T(i_T, j_T+1) - T(i_T, j_T-1)) / (2 * dy(j_cv));
    end
end

%  k at  nodes
k_n = zeros(Nx_n, Ny_n);
for i = 1:Nx_n
    for j = 1:Ny_n
        k_n(i, j) = k(X_n(i), Y_n(j));
   
    end
end

qx = - (k_n .* dTdx); 
qy = - (k_n .* dTdy);


figure(3);
contourf(X, Y, T', 25);hold on;
colorbar;
quiver(X, Y, qx', qy','Color','red'); 
axis equal tight;
title('Heat Flux Vector Plot ');
xlabel('x');
ylabel('y');

% Print some useful info for rough comp
T_max = max(T(:));
T_min =  min(T(:));

% index of centreline in x-direction and y dirn ---
x_mid = max(X_n(:)) / 2   ;                     
[~, i_cen] = min(abs(X_n(:,1) - x_mid));     
y_mid = max(Y_n(:)) / 2   ;                     
[~, j_cen] = min(abs(Y_n(:,1) - y_mid));
T_centl = T(i_cen,   :);   
T_avg_centl = mean(T_centl );

fprintf(' Average temperatu re along vertical centreline = %.4f\n'  , T_avg_centl);


fprintf('T_max = %.6f, T_min = %.6f\n', T_max, T_min);
fprintf("----------------------------------------------------------------------\n")
% Optionally show a few cross-sections
figure(4);
plot(T(i_cen, :) ,Y_n,  '-o','DisplayName',sprintf("%d X %d" , Nx,Ny)); hold on; % centerline vs y 
ylabel('y');
xlabel('T (centerline y)');
title('Centerline Temperature vs y');
legend show;

figure(5);
plot(X_n, T(:,j_cen) , '-o','DisplayName',sprintf("%d X %d ", Nx ,Ny) );hold on; % centerline vs x
xlabel('x');
ylabel('T (centerline y)');
title('Centerline Temperature vs x ');
legend show;




figure(6);
[X1, Y1] = meshgrid(X_n, Y_n);
[x1 , y1 ] = meshgrid( gx , gy ) ;

plot(x1, gy, color = 'blue') ; hold on;
plot(gx, y1' , color = 'blue') ;
xlabel('Lx')
ylabel("Ly")
%Plotting the node points
plot(X1, Y_n,'r.');



%save(sprintf('T %d X %d, grx = %0.2f .mat',  Nx,Ny,rx),"T","X_n","Y_n","qx",'qy');