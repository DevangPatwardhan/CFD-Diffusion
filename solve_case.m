function solve_case(caseID, L, H, meshList, tol, maxIter, omega, outdir)
% solve_case: sets up BCs, runs meshes and plots results
if nargin<8, outdir = pwd; end


% get case definitions
[BC, S_fun, k_fun, caseName] = case_definitions(caseID, L, H);


% arrays to keep solutions for mesh independence
solutions = cell(numel(meshList),1);
xcells = cell(numel(meshList),1);
ycells = cell(numel(meshList),1);
resids = cell(numel(meshList),1);
iters = zeros(numel(meshList),1);


for m = 1:numel(meshList)
nx = meshList(m); ny = meshList(m);
fprintf('\n-- Mesh %dx%d --\n', nx, ny);
% choose stretch factors (example: cluster near x=0 if BC has strong gradients)
stretchX = 0; stretchY = 0;
[xc, yc] = grid_gen(L, H, nx, ny, stretchX, stretchY);
% initial guess: linear in x between Dirichlet if available
Tinit = zeros(ny, nx);
Tinit(:) = mean([0,100]);
[T, residuals, nit] = fvm_solver(Tinit, xc, yc, k_fun, S_fun, BC, tol, maxIter, omega);
solutions{m} = T; xcells{m}=xc; ycells{m}=yc; resids{m}=residuals; iters(m)=nit;


% plotting and saving
fname = fullfile(outdir, sprintf('case%d_mesh%dx%d.mat', caseID, nx, ny));
save(fname, 'T','xc','yc','residuals','nit','BC','caseID');


% plot temperature and flux
figure('Visible','off');
utils_plot.contour_with_flux(xc, yc, T, k_fun, sprintf('Case %d (%s) - T, %dx%d', caseID, caseName, nx, ny));
saveas(gcf, fullfile(outdir, sprintf('case%d_T_mesh%dx%d.png',caseID,nx,ny)));
close;


figure('Visible','off'); semilogy(residuals); xlabel('Iteration'); ylabel('Residual (max |dT|)'); title(sprintf('Residual history - Case %d mesh %dx%d',caseID, nx,ny)); grid on;
saveas(gcf, fullfile(outdir, sprintf('case%d_res_mesh%dx%d.png',caseID,nx,ny)));
close;
end


% mesh independence: centerline y = H/2
figure('Visible','off'); hold on;
xq = linspace(0,L,200);
for m = 1:numel(meshList)
Tline = interp2(xcells{m}, ycells{m}, solutions{m}, xq, (H/2)*ones(size(xq)));
plot(xq, Tline, 'DisplayName', sprintf('%dx%d', meshList(m), meshList(m)));
end
xlabel('x'); ylabel('T'); title(sprintf('Mesh independence - Case %d: %s', caseID, caseName)); legend show; grid on;
saveas(gcf, fullfile(outdir, sprintf('case%d_mesh_independence.png',caseID)));
close;


% save summary
summary.caseID = caseID; summary.caseName = caseName; summary.meshList = meshList; summary.iters = iters;
save(fullfile(outdir, sprintf('case%d_summary.mat',caseID)),'summary');


fprintf('\nFinished case %d. Saved results to %s\n', caseID, outdir);
end