function [xc, yc] = grid_gen(L, H, nx, ny, stretchX, stretchY)
% grid_gen: returns cell-centered coordinates for non-uniform grid
if nargin<6, stretchY = 0; end
if nargin<5, stretchX = 0; end


xc = gen1d(L, nx, stretchX);
yc = gen1d(H, ny, stretchY);


function x = gen1d(len, n, stretch)
if stretch == 0
dx = len/n; x = dx/2 : dx : len-dx/2;
else
beta = stretch; s = (0.5:1:n-0.5)/n; xi = (exp(beta*s)-1)/(exp(beta)-1); x = xi * len;
end
end
end