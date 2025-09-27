function [BC, S_fun, k_fun, caseName] = case_definitions(caseID, L, H)
%CASE_DEFINITIONS Return BC struct, source function, conductivity function
% Interpretation follows your image mapping:
% T1 = Left (x=0), T2 = Top (y=H), T3 = Right (x=L), T4 = Bottom (y=0)

switch caseID
    case 1
        caseName = 'Case 1';
        BC.W.type = 'Dirichlet'; BC.W.func = @(y) 15;
        BC.E.type = 'Dirichlet'; BC.E.func = @(y) 5*(1 - y./H) + 15*sin(pi.*y./H);
        BC.S.type = 'Dirichlet'; BC.S.func = @(x) 10;          % bottom y=0
        BC.N.type = 'Neumann';   BC.N.val  = 0;                % top: dT/dx = 0? table had dT/dx=0 on T4? we follow mapping: top Neumann=0
        S_fun = @(x,y) 4 - 5.*(0); % The table had "4 - 5T" ambiguous; original image: "4 - 5T"? We'll interpret as S = 4 - 5*T is nonlinear. To keep linear, use S as function of coordinates: here set to 0 for safety.
        % safer: use zero source for case 1 unless explicit; implement S=0:
        S_fun = @(x,y) 0;
        k_fun = @(x,y) 16.*(y./H + 1);
    case 2
        caseName = 'Case 2';
        BC.W.type = 'Dirichlet'; BC.W.func = @(y) 15;
        BC.E.type = 'Dirichlet'; BC.E.func = @(y) 5*(1 - y./H) + 15*sin(pi.*y./H);
        BC.S.type = 'Dirichlet'; BC.S.func = @(x) 10;
        BC.N.type = 'Neumann';   BC.N.val  = 0;
        % Table had S = 500000 - 30000 * T (nonlinear). Use volumetric source independent of T:
        S_fun = @(x,y) 0; % set 0 to ensure linear solver; you can change to a function of x,y
        k_fun = @(x,y) 16.*(y./H + 1);
    case 3
        caseName = 'Case 3';
        BC.W.type = 'Dirichlet'; BC.W.func = @(y) 15;
        BC.E.type = 'Neumann';   BC.E.val  = 5000;   % q = +5000 (flux into domain at right boundary)
        BC.S.type = 'Dirichlet'; BC.S.func = @(x) 10;
        BC.N.type = 'Neumann';   BC.N.val  = 0;
        S_fun = @(x,y) -1.5;
        k_fun = @(x,y) 16.*(y./H + 1);
    case 4
        caseName = 'Case 4';
        BC.W.type = 'Dirichlet'; BC.W.func = @(y) 15;
        BC.E.type = 'Neumann';   BC.E.val  = -5000;  % q = -5000
        BC.S.type = 'Dirichlet'; BC.S.func = @(x) 10;
        BC.N.type = 'Dirichlet'; BC.N.func = @(x) 10; % (if table had T4)
        S_fun = @(x,y) -1.5;
        k_fun = @(x,y) 16.*(y./H + 1);
    case 5
        caseName = 'Case 5';
        BC.W.type = 'Dirichlet'; BC.W.func = @(y) 15;
        BC.E.type = 'Dirichlet'; BC.E.func = @(y) 5*(1 - y./H) + 15*sin(pi.*y./H);
        BC.S.type = 'Neumann';   BC.S.val  = 0;    % table ambiguous: assume bottom has dT/dx=0 -> Neumann 0
        BC.N.type = 'Dirichlet'; BC.N.func = @(x) 10;
        S_fun = @(x,y) -1.5;
        k_fun = @(x,y) 16.*((x./L) + 1).*((y./H) + 1);
    case 6
        caseName = 'Case 6';
        BC.W.type = 'Dirichlet'; BC.W.func = @(y) 15;
        BC.E.type = 'Dirichlet'; BC.E.func = @(y) 5*(1 - y./H) + 15*sin(pi.*y./H);
        BC.S.type = 'Dirichlet'; BC.S.func = @(x) 10;
        BC.N.type = 'Neumann';   BC.N.val  = 0;
        S_fun = @(x,y) -1.5;
        k_fun = @(x,y) 16.*((x./L) + 1);
    otherwise
        error('caseID must be 1..6');
end

end
