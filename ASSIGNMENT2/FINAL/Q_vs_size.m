%  NOTE : CHATGPT Was used for easy plotting of these quantities :)
mesh_sizes = [20 30 40 60 80];

% Temperatures
T_avg = [11.4415 11.4637 11.3597 11.2926 11.3799];
T_max = [15 15 15 15 15];
T_min = [-17.537172 -17.348445 -17.268244 -17.196071 -17.162387];

% Solver performance
iterations = [1245 3044 6154 20363 61249];  % number of iterations
time_taken = [0.093867 0.412913 1.445961 10.598615 73.177766]; % in seconds

% --- Plot Average Temperature ---
figure;
plot(mesh_sizes, T_avg, '-o','LineWidth',1.5,'MarkerSize',8);
xlabel('Number of cells (N_x = N_y)');
ylabel('Average Temperature');
title('Average Temperature vs Mesh Size');
grid on;

% --- Plot Maximum Temperature ---
figure;
plot(mesh_sizes, T_max, '-s','LineWidth',1.5,'MarkerSize',8);
xlabel('Number of cells (N_x = N_y)');
ylabel('Maximum Temperature');
title('Maximum Temperature vs Mesh Size');
grid on;

% --- Plot Minimum Temperature ---
figure;
plot(mesh_sizes, T_min, '-^','LineWidth',1.5,'MarkerSize',8);
xlabel('Number of cells (N_x = N_y)');
ylabel('Minimum Temperature');
title('Minimum Temperature vs Mesh Size');
grid on;

% --- Plot Iterations ---
figure;
plot(mesh_sizes, iterations, '-d','LineWidth',1.5,'MarkerSize',8);
xlabel('Number of cells (N_x = N_y)');
ylabel('Iterations to Converge');
title('Iterations vs Mesh Size');
grid on;

% --- Plot Time Taken ---
figure;
plot(mesh_sizes, time_taken, '-p','LineWidth',1.5,'MarkerSize',8);
xlabel('Number of cells (N_x = N_y)');
ylabel('Time Taken (s)');
title('Time Taken vs Mesh Size');
grid on;
