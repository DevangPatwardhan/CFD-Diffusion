function run_all()
%RUN_ALL top-level runner - loads params and executes the chosen case over meshList
    P = params();
    if ~exist(P.outdir,'dir')
        mkdir(P.outdir);
    end

    [BC, S_fun, k_fun, caseName] = case_definitions(P.caseID, P.L, P.H);

    for m = 1:numel(P.meshList)
        meshCells = P.meshList(m); % cell count
        Nx_nodes = meshCells + 1; Ny_nodes = meshCells + 1;
        fprintf('\n-- Running Case %d (%s) mesh %dx%d --\n', P.caseID, caseName, meshCells, meshCells);

        % build mesh
        [xc, yc, x_nodes, y_nodes, dx_cell, dy_cell, nx, ny] = generate_mesh(Nx_nodes, Ny_nodes, P.rx, P.ry, P.L, P.H);

        % initial guess
        Tinit = ones(ny, nx) * 15;

        % solve
        [T, residuals, nit] = fvm_solver(Tinit, xc, yc, k_fun, S_fun, BC, P.tol, P.maxIter, P.omega);

        % save data
        fname = fullfile(P.outdir, sprintf('case%d_mesh%dx%d.mat', P.caseID, meshCells, meshCells));
        save(fname, 'T','xc','yc','x_nodes','y_nodes','residuals','nit','P');

        % plots
        figtitle = sprintf('Case %d (%s) mesh %dx%d', P.caseID, caseName, meshCells, meshCells);
        outpng = fullfile(P.outdir, sprintf('case%d_T_mesh%dx%d.png',P.caseID,meshCells,meshCells));
        utils_plot_contour_with_flux(xc, yc, T, k_fun, figtitle, outpng);

        % residual plot
        figure('Visible','off'); semilogy(1:numel(residuals), residuals, '-o'); grid on;
        xlabel('Iteration'); ylabel('Residual (max |dT|)'); title(['Residual history - ' figtitle]);
        saveas(gcf, fullfile(P.outdir, sprintf('case%d_res_mesh%dx%d.png',P.caseID,meshCells,meshCells)));
        close;

    end

    % mesh independence centerline plot (y = H/2)
    figure;
    hold on;
    xq = linspace(0, P.L, 200);
    for m = 1:numel(P.meshList)
        meshCells = P.meshList(m);
        load(fullfile(P.outdir, sprintf('case%d_mesh%dx%d.mat', P.caseID, meshCells, meshCells)), 'T','xc','yc');
        Tline = interp2(xc, yc, T, xq, (P.H/2)*ones(size(xq)));
        plot(xq, Tline, 'DisplayName', sprintf('%dx%d', meshCells, meshCells));
    end
    xlabel('x'); ylabel('T'); title(sprintf('Mesh independence - Case %d',P.caseID)); legend show; grid on;
    saveas(gcf, fullfile(P.outdir, sprintf('case%d_mesh_independence.png',P.caseID)));
    close;

    fprintf('\nAll done. Results in folder: %s\n', P.outdir);
end
