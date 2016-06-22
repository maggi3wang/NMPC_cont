function u_aux = aux_control(rho,X,X_dot,K,w,B,W_fnc,m)

%X: geodesic
%X_dot: geodesic direction
%K: number of points along geodesic
%w: integration quadrature weights

% rho = @(x) 2+1.298477064463527*x(1) +...
%             +9.936724981692999*(x(1)^2);
% rho = @(x) 2+2.821996084008429*x(1) +...
%             +4.422187495932881e+03*(x(1)^2);

u_aux = zeros(m,1);

for k = 1:K+1
    x_k = X(:,k);
    W = W_fnc(x_k);
    x_dot_k = X_dot(:,k);
    u_aux = u_aux + (1/2)*w(k)*(-1/2)*rho(x_k)*B'*(W\x_dot_k);
end
    



end