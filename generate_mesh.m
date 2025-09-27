function [xc, yc, x_nodes, y_nodes, dx_cell, dy_cell, nx, ny] = generate_mesh(Nx, Ny, rx, ry, L, H)
%GENERATE_MESH Generate stretched nodes and cell-centers for FVM.
% Inputs:
%   Nx,Ny = number of nodes in x,y (including boundary nodes)
%   rx,ry = stretching factors in x,y (your formula)
%   L,H   = domain sizes
% Outputs:
%   xc,yc        = cell-center coordinates (vectors)
%   x_nodes,y_nodes = node coords (including boundaries)
%   dx_cell,dy_cell = cell widths arrays (length nx,ny)
%   nx, ny       = number of cells in x,y (length(xc),length(yc))

    % generate nodes (0..L and 0..H) using user's algorithm
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

    % cell centers are midpoints between adjacent nodes
    xc = 0.5 * (x_nodes(1:end-1) + x_nodes(2:end));
    yc = 0.5 * (y_nodes(1:end-1) + y_nodes(2:end));

    nx = numel(xc);
    ny = numel(yc);

    % compute local cell widths (distance between adjacent faces)
    dx_face = diff(x_nodes); % face widths
    dy_face = diff(y_nodes);

    % cell widths (cell i sits between face i and i+1) -> dx_cell = dx_face
    dx_cell = dx_face; % length nx
    dy_cell = dy_face; % length ny

end
