% I have used AI for repetitive coding here

% Load 10x10
D10 = load("T_ eq 10 X 10.mat");
T10  = D10.T;
Xn10 = D10.X_n;
Yn10 = D10.Y_n;
qx10 = D10.qx;
qy10 = D10.qy;

% Load 20x20
D20 = load("T_ eq 20 X 20.mat");
T20  = D20.T;
Xn20 = D20.X_n;
Yn20 = D20.Y_n;
qx20 = D20.qx;
qy20 = D20.qy;

% Load 30x30
D30 = load("T_ eq 30 X 30.mat");
T30  = D30.T;
Xn30 = D30.X_n;
Yn30 = D30.Y_n;
qx30 = D30.qx;
qy30 = D30.qy;

% Load 40x40
D40 = load("T_ eq 40 X 40.mat");
T40  = D40.T;
Xn40 = D40.X_n;
Yn40 = D40.Y_n;
qx40 = D40.qx;
qy40 = D40.qy;

% Load 80x80
D80 = load("T_ 80 X 80.mat");
T80  = D80.T;
Xn80 = D80.X_n;
Yn80 = D80.Y_n;
qx80 = D80.qx;
qy80 = D80.qy;


% Load all variables
data10 = load("T_ eq 10 X 10.mat"); 
T10 = data10.T; Xn10 = data10.X_n; Yn10 = data10.Y_n;

data20 = load("T_ eq 20 X 20.mat"); 
T20 = data20.T; Xn20 = data20.X_n; Yn20 = data20.Y_n;

data30 = load("T_ eq 30 X 30.mat"); 
T30 = data30.T; Xn30 = data30.X_n; Yn30 = data30.Y_n;

data40 = load("T_ eq 40 X 40.mat"); 
T40 = data40.T; Xn40 = data40.X_n; Yn40 = data40.Y_n;

data80 = load("T_ 80 X 80.mat"); 
T80 = data80.T; Xn80 = data80.X_n; Yn80 = data80.Y_n;
% Plot all centerlines
figure; hold on;

Ny10 = length(Yn10);
plot(Xn10, T10(:, round(Ny10/2)), '-o', 'DisplayName','10x10');

Ny20 = length(Yn20);
plot(Xn20, T20(:, round(Ny20/2)), '-s', 'DisplayName','20x20');

Ny30 = length(Yn30);
plot(Xn30, T30(:, round(Ny30/2)), '-^', 'DisplayName','30x30');

Ny40 = length(Yn40);
plot(Xn40, T40(:, round(Ny40/2)), '-d', 'DisplayName','40x40');

Ny80 = length(Yn80);
plot(Xn80, T80(:, round(Ny80/2)), '-d', 'DisplayName','80x80');

xlabel('x');
ylabel('T (centerline y)');
title('Centerline Temperature vs x');
legend show;
grid on;

% Domain midpoints (adjust if L,H are not 1.0)
x_mid = max(Xn10)/2;
y_mid = max(Yn10)/2;

% --- For 10x10 mesh ---
[~, ix10] = min(abs(Xn10 - x_mid));
[~, iy10] = min(abs(Yn10 - y_mid));
Tc10 = T10(ix10, iy10);

% --- For 20x20 mesh ---
x_mid = max(Xn20)/2;  % (ensure consistent domain size)
y_mid = max(Yn20)/2;
[~, ix20] = min(abs(Xn20 - x_mid));
[~, iy20] = min(abs(Yn20 - y_mid));
Tc20 = T20(ix20, iy20);

% --- For 30x30 mesh ---
x_mid = max(Xn30)/2;
y_mid = max(Yn30)/2;
[~, ix30] = min(abs(Xn30 - x_mid));
[~, iy30] = min(abs(Yn30 - y_mid));
Tc30 = T30(ix30, iy30);

% --- For 40x40 mesh ---
x_mid = max(Xn40)/2;
y_mid = max(Yn40)/2;
[~, ix40] = min(abs(Xn40 - x_mid));
[~, iy40] = min(abs(Yn40 - y_mid));
Tc40 = T40(ix40, iy40);

% --- For 80x80 mesh ---
x_mid = max(Xn80)/2;
y_mid = max(Yn80)/2;
[~, ix80] = min(abs(Xn80 - x_mid));
[~, iy80] = min(abs(Yn80 - y_mid));
Tc80 = T80(ix80, iy80);

% Collect results including 80x80
mesh_sizes = [10 20 30 40 80];
Tmid = [Tc10 Tc20 Tc30 Tc40 Tc80];

% Plot
figure;
plot(mesh_sizes, Tmid, '-o','LineWidth',1.5,'MarkerSize',8);
xlabel('Number of cells (N_x = N_y)');
ylabel('Temperature at geometric midpoint');
title('Midpoint Temperature vs Mesh Size (Stretched Grid)');
grid on;


% Plot
figure;
plot(mesh_sizes, Tmid, '-o','LineWidth',1.5,'MarkerSize',8);
xlabel('Number of cells (N_x = N_y)');
ylabel('Temperature at geometric midpoint');
title('Midpoint Temperature vs Mesh Size (Stretched Grid)');
grid on;




