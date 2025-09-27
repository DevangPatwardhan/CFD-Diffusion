function utils_plot_contour_with_flux(xc, yc, T, k_fun, titleStr, outpath)
% Simple helper: contour + quiver overlay and optional save
    [Xc, Yc] = meshgrid(xc, yc);
    figure('Visible','off');
    contourf(Xc, Yc, T, 30, 'LineStyle','none'); colorbar;
    hold on;
    [qx, qy] = flux_compute(T, xc, yc, k_fun);
    skip = max(1, round(numel(xc)/20));
    quiver(Xc(1:skip:end,1:skip:end), Yc(1:skip:end,1:skip:end), qx(1:skip:end,1:skip:end), qy(1:skip:end,1:skip:end));
    axis equal tight; xlabel('x'); ylabel('y'); title(titleStr);
    if nargin>=6 && ~isempty(outpath)
        saveas(gcf, outpath);
    end
    close;
end
