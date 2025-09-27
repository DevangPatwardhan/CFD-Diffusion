function [qx, qy] = flux_compute(T, xc, yc, k_fun)
% Compute q = -k grad T at cell centers (finite differences)
    nx = numel(xc); ny = numel(yc);
    qx = zeros(ny, nx); qy = zeros(ny, nx);

    for j=1:ny
        for i=1:nx
            % dT/dx
            if i==1
                dx = xc(2)-xc(1);
                dTdx = (T(j,2)-T(j,1))/dx;
            elseif i==nx
                dx = xc(end)-xc(end-1);
                dTdx = (T(j,end)-T(j,end-1))/dx;
            else
                dx = xc(i+1)-xc(i-1);
                dTdx = (T(j,i+1)-T(j,i-1))/dx;
            end
            % dT/dy
            if j==1
                dy = yc(2)-yc(1);
                dTdy = (T(2,i)-T(1,i))/dy;
            elseif j==ny
                dy = yc(end)-yc(end-1);
                dTdy = (T(end,i)-T(end-1,i))/dy;
            else
                dy = yc(j+1)-yc(j-1);
                dTdy = (T(j+1,i)-T(j-1,i))/dy;
            end
            kP = k_fun(xc(i), yc(j));
            qx(j,i) = -kP * dTdx;
            qy(j,i) = -kP * dTdy;
        end
    end
end
