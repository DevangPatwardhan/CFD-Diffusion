function P = params()
%PARAMS Define all parameters for the diffusion FVM solver (single central file).
    % Domain
    P.L = 1.0;     % Length in x
    P.H = 0.5;     % Height in y

    % Mesh settings (number of NODES, so cells = Nx-1, Ny-1)
    P.Nx = 11;     % Number of nodes in x (including boundaries)
    P.Ny = 11;     % Number of nodes in y (including boundaries)
    P.rx = 1.2;    % Stretching factor in x (your script style)
    P.ry = 1.2;    % Stretching factor in y

    % Solver controls
    P.tol     = 1e-6;   % Convergence tolerance
    P.maxIter = 10000;  % Maximum iterations
    P.omega   = 1.0;    % Relaxation factor (1 = GS, >1 = SOR)

    % Case selection (1..6 per your table image)
    P.caseID = 4;

    % Mesh independence study (list of cell counts, e.g. nodes-1)
    P.meshList = [10, 20, 40];  % will use Nx = mesh+1 nodes internally

    % Output folder
    P.outdir = 'results';

end
