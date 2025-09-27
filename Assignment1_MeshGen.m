function [xc, yc, dx, dy, nx, ny] = generate_mesh(Nx, Ny, rx, ry, L, H)
%GENERATE_MESH  Generate stretched grid and return cell centers.
%
% Inputs:
%   Nx, Ny  - number of nodes in x, y (including boundaries)
%   rx, ry  - stretching factors
%   L, H    - domain lengths
%
% Outputs:
%   xc, yc  - vectors of cell center coords
%   dx, dy  - average cell widths
%   nx, ny  - number of cells in x, y

    % Generate stretched grid nodes in [0,L] and [0,H]
    x_nodes = zeros(1, Nx);
    x_nodes(1) = 0;
    for i = 2:Nx
        x_nodes(i) = x_nodes(i-1) + (1 - x_nodes(i-1)) * (1 - 1/rx);
    end
    x_nodes = x_nodes / x_nodes(end) * L;

    y_nodes = zeros(1, Ny);
    y_nodes(1) = 0;
    for i = 2:Ny
        y_nodes(i) = y_nodes(i-1) + (1 - y_nodes(i-1)) * (1 - 1/ry);
    end
    y_nodes = y_nodes / y_nodes(end) * H;

    % Cell centers (midpoints between nodes)
    xc = 0.5 * (x_nodes(1:end-1) + x_nodes(2:end));
    yc = 0.5 * (y_nodes(1:end-1) + y_nodes(2:end));

    % Number of cells
    nx = length(xc);
    ny = length(yc);

    % Approx cell sizes (uniform avg)
    dx = mean(diff(x_nodes));
    dy = mean(diff(y_nodes));

end