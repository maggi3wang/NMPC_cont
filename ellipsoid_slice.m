function [alpha]= ellipsoid_slice(M, RPI, p_coords,h_coords,h_vals,h_constr,center,N,c)

%M: metric
%RPI: bound
%p_coords: plane coordinates
%h_coords: hidden coordinates
%h_vals: values of hidden coordinates

n = size(M,1);

gamma = linspace(0,2*pi,N);
gamma = gamma';

alpha = zeros(N,1);

%%
if (h_constr)
    for i = 1:N
        cvx_begin quiet
        variables a v(n)
        maximize (a)
        subject to
        x = a*cos(gamma(i));
        y = a*sin(gamma(i));
        a >= 0;
        v(p_coords) == center(p_coords) + [x;y];
        v(h_coords) == h_vals;
        (v-center)'*M*(v-center) <= RPI;
        cvx_end
        if strcmp(cvx_status,'Solved') == 1
            alpha(i) = a;
        end
    end
else
    for i = 1:N
        cvx_begin quiet
        variables a v(n)
        maximize (a)
        subject to
        x = a*cos(gamma(i));
        y = a*sin(gamma(i));
        a >= 0;
        v(p_coords) == center(p_coords) + [x;y];
        (v-center)'*M*(v-center) <= RPI;
        cvx_end
        if strcmp(cvx_status,'Solved') == 1
            alpha(i) = a;
        end
    end
end
%%
ellipse_x = center(p_coords(1)) + alpha.*cos(gamma);
ellipse_y = center(p_coords(2)) + alpha.*sin(gamma);

plot(ellipse_x,ellipse_y,c,'linewidth',2);
hold on
patch(ellipse_x,ellipse_y,c,'facealpha',0.1);
grid on
axis equal
end

