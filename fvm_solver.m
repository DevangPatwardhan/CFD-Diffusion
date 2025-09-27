function [T, residuals, nit, info] = fvm_solver_newton(Tinit, xc, yc, k_fun, S_fun, BC, tol, maxIter, omega)
%FVM_SOLVER_NEWTON  Nonlinear FVM solver using Newton-Raphson + sparse linear solves.
%
% [T, residuals, nit, info] = fvm_solver_newton(Tinit, xc, yc, k_fun, S_fun, BC, tol, maxIter, omega)
%
% Inputs:
%   Tinit : ny x nx initial guess
%   xc,yc : cell-center coordinate vectors
%   k_fun : k(x,y) conductivity handle
%   S_fun : source handle. Preferably [S,dSdT] = S_fun(x,y,T). Fallbacks:
%           S = S_fun(x,y,T) or S = S_fun(x,y).
%   BC    : boundary struct with fields W,E,S,N. For Dirichlet use .type='Dirichlet' & .func=@(coord)Tbc
%           For Neumann use .type='Neumann' & .val = flux (positive = into domain)
%   tol,maxIter,omega : Newton tolerance, maximum Newton iterations, step damping (0<omega<=1)
%
% Outputs:
%   T         : ny x nx solution
%   residuals : vector of ||delta||_inf per Newton iteration
%   nit       : number of Newton iterations performed
%   info      : struct with .linSolves and .finalResidual

if nargin < 9, omega = 1.0; end
if nargin < 8, maxIter = 200; end
if nargin < 7, tol = 1e-6; end

% grid sizes
nx = numel(xc); ny = numel(yc); nc = nx * ny;

T = Tinit;
residuals = [];
linSolves = 0;

% Build faces from centers (extrapolate)
x_face = zeros(1, nx+1);
if nx>1
    x_face(2:nx) = 0.5*(xc(1:end-1) + xc(2:end));
    x_face(1) = xc(1) - (x_face(2) - xc(1));
    x_face(end) = xc(end) + (xc(end) - x_face(end-1));
else
    x_face = [0, xc(1)*2];
end
y_face = zeros(1, ny+1);
if ny>1
    y_face(2:ny) = 0.5*(yc(1:end-1) + yc(2:end));
    y_face(1) = yc(1) - (y_face(2) - yc(1));
    y_face(end) = yc(end) + (yc(end) - y_face(end-1));
else
    y_face = [0, yc(1)*2];
end

% cell widths between faces
dx_cell = diff(x_face); % length nx
dy_cell = diff(y_face); % length ny

% cell areas (ny x nx): meshgrid gives rows = dy (ny) and cols = dx (nx)
[DX, DY] = meshgrid(dx_cell, dy_cell); % size ny x nx
cellArea = DX .* DY;

% helper index mapping p = i + (j-1)*nx  (i: 1..nx, j:1..ny)
idx = @(i,j) i + (j-1)*nx;

% Newton loop
for it = 1:maxIter
    % Preallocate sparse assembly arrays (max 5 entries per cell)
    maxEntries = 5 * nc;
    I = zeros(maxEntries,1); J = zeros(maxEntries,1); V = zeros(maxEntries,1); ctr = 0;
    R = zeros(nc,1); % residual vector
    
    for j = 1:ny
        for i = 1:nx
            p = idx(i,j);
            xp = xc(i); yp = yc(j);
            
            % face areas for cell (2D)
            Aw = dy_cell(j); Ae = dy_cell(j);
            As = dx_cell(i); An = dx_cell(i);
            
            % distances center-to-center or center-to-boundary-face
            if i>1, dPW = xc(i) - xc(i-1); else dPW = xc(i) - x_face(1); end
            if i<nx, dPE = xc(i+1) - xc(i); else dPE = x_face(end) - xc(i); end
            if j>1, dPS = yc(j) - yc(j-1); else dPS = yc(j) - y_face(1); end
            if j<ny, dPN = yc(j+1) - yc(j); else dPN = y_face(end) - yc(j); end
            
            % evaluate k at P and neighbors (where needed)
            kP = k_fun(xp, yp);
            
            % initialize conductances and source terms
            aW = 0; aE = 0; aS = 0; aN = 0; Su = 0; Sp = 0;
            
            % WEST face
            if i>1
                kW = k_fun(xc(i-1), yp);
                kf = 2/(1/kP + 1/kW);
                aW = kf * Aw / dPW;
            else
                if strcmpi(BC.W.type, 'Dirichlet')
                    Tb = BC.W.func(yp);
                    dface = xc(i) - x_face(1);
                    kf = kP;
                    aW = kf * Aw / dface;
                    Su = Su + aW * Tb;  % move Dirichlet to source
                    aW = 0;
                else
                    qn = BC.W.val; % positive into domain
                    Su = Su + qn * Aw;
                end
            end
            
            % EAST face
            if i<nx
                kE = k_fun(xc(i+1), yp);
                kf = 2/(1/kP + 1/kE);
                aE = kf * Ae / dPE;
            else
                if strcmpi(BC.E.type, 'Dirichlet')
                    Tb = BC.E.func(yp);
                    dface = x_face(end) - xc(i);
                    kf = kP;
                    aE = kf * Ae / dface;
                    Su = Su + aE * Tb;
                    aE = 0;
                else
                    qn = BC.E.val; % positive into domain
                    Su = Su + qn * Ae;
                end
            end
            
            % SOUTH face
            if j>1
                kS = k_fun(xp, yc(j-1));
                kf = 2/(1/kP + 1/kS);
                aS = kf * As / dPS;
            else
                if strcmpi(BC.S.type, 'Dirichlet')
                    Tb = BC.S.func(xp);
                    dface = yc(j) - y_face(1);
                    kf = kP;
                    aS = kf * As / dface;
                    Su = Su + aS * Tb;
                    aS = 0;
                else
                    qn = BC.S.val; % positive into domain
                    Su = Su + qn * As;
                end
            end
            
            % NORTH face
            if j<ny
                kN = k_fun(xp, yc(j+1));
                kf = 2/(1/kP + 1/kN);
                aN = kf * An / dPN;
            else
                if strcmpi(BC.N.type, 'Dirichlet')
                    Tb = BC.N.func(xp);
                    dface = y_face(end) - yc(j);
                    kf = kP;
                    aN = kf * An / dface;
                    Su = Su + aN * Tb;
                    aN = 0;
                else
                    qn = BC.N.val; % positive into domain
                    Su = Su + qn * An;
                end
            end
            
            % Evaluate source S and derivative dS/dT at (xp,yp,Tp)
            Tp = T(j,i);
            % Try preferred signature [S,dSdT] = S_fun(x,y,T)
            dSdT = 0;
            try
                [Sval, dSdT] = S_fun(xp, yp, Tp);
            catch ME1
                % If not supported, try S_fun(x,y,T) returning S only
                try
                    Sval = S_fun(xp, yp, Tp);
                    % numeric derivative
                    epsT = max(1e-8, 1e-6 * max(1, abs(Tp)));
                    Splus = S_fun(xp, yp, Tp + epsT);
                    dSdT = (Splus - Sval) / epsT;
                catch ME2
                    % last fallback: S_fun(x,y) (T-independent)
                    Sval = S_fun(xp, yp);
                    dSdT = 0;
                end
            end
            
            % integrated source and Jacobian contribution
            area = cellArea(j,i);
            Sc = Sval * area;
            Sp = dSdT * area; % Jacobian addition (positive if dS/dT>0)
            
            % diagonal coefficient (include Sp with PLUS sign)
            aP = aW + aE + aS + aN + Sp;
            
            % neighbor contribution sum (for residual)
            neighSum = 0;
            if i>1 && aW>0, neighSum = neighSum + aW * T(j, i-1); end
            if i<nx && aE>0, neighSum = neighSum + aE * T(j, i+1); end
            if j>1 && aS>0, neighSum = neighSum + aS * T(j-1, i); end
            if j<ny && aN>0, neighSum = neighSum + aN * T(j+1, i); end
            
            % Residual for cell p (to satisfy A*delta = -R)
            R(p) = aP * Tp - neighSum - Su - Sc;
            
            % Assemble sparse entries: diag + neighbors (-aW etc.)
            ctr = ctr + 1; I(ctr) = p; J(ctr) = p; V(ctr) = aP;
            if i>1 && aW>0
                ctr = ctr + 1; I(ctr) = p; J(ctr) = idx(i-1, j); V(ctr) = -aW;
            end
            if i<nx && aE>0
                ctr = ctr + 1; I(ctr) = p; J(ctr) = idx(i+1, j); V(ctr) = -aE;
            end
            if j>1 && aS>0
                ctr = ctr + 1; I(ctr) = p; J(ctr) = idx(i, j-1); V(ctr) = -aS;
            end
            if j<ny && aN>0
                ctr = ctr + 1; I(ctr) = p; J(ctr) = idx(i, j+1); V(ctr) = -aN;
            end
        end
    end
    
    % truncate assembly arrays
    I = I(1:ctr); J = J(1:ctr); V = V(1:ctr);
    A = sparse(I, J, V, nc, nc);
    
    % Solve linear system A * delta = -R
    linSolves = linSolves + 1;
    % guard against singular or ill-conditioned A
    try
        delta = A \ (-R);
    catch ME
        warning('fvm_solver_newton:linsolveFail', 'Linear solve failed at Newton it %d: %s', it, ME.message);
        nit = it;
        info.linSolves = linSolves;
        info.finalResidual = norm(R, Inf);
        return;
    end
    
    % damping / step scale
    delta = omega * delta;
    
    % update T (reshape delta into nx-by-ny, then transpose -> ny-by-nx)
    % delta vector ordered by p = i + (j-1)*nx, so reshape into [nx, ny] then transpose
    dT = reshape(delta, [nx, ny]).';
    T = T + dT;
    
    delta_norm = norm(delta, Inf);
    residuals(end+1) = delta_norm;
    
    % check convergence
    if delta_norm < tol
        nit = it;
        info.linSolves = linSolves;
        info.finalResidual = delta_norm;
        return;
    end
    
    % safety checks for divergence
    if any(isnan(T(:))) || any(isinf(T(:))) || delta_norm > 1e12
        warning('fvm_solver_newton:diverge', 'Newton diverged at it=%d (||delta||=%g).', it, delta_norm);
        nit = it;
        info.linSolves = linSolves;
        info.finalResidual = delta_norm;
        return;
    end
end

% did not converge within maxIter
nit = maxIter;
info.linSolves = linSolves;
if isempty(residuals), info.finalResidual = NaN; else info.finalResidual = residuals(end); end
warning('fvm_solver_newton:noConverge', 'Newton did not converge in %d iterations (last ||delta||=%g)', maxIter, info.finalResidual);

end
