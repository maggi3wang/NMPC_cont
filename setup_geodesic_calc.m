function [geo_Prob,K_e,w_e,T_e,T_dot_e,Aeq] = setup_geodesic_calc(n,N,W,dW)

K = N+4;
K_e = K+2;

%Obtain Chebyschev Pseudospectral Numerics

%CGL points and quadrature weights
[t,w] = clencurt(K);

[t_e,w_e] = clencurt(K_e);

% Chebyshev polynomial method
[phi_start, phi_dot_start] = compute_cheby(0,N,-1);
[phi_end, phi_dot_end] = compute_cheby(0,N,1);
A_start = kron(eye(n),phi_start');
A_end = kron(eye(n),phi_end');

Aeq = [A_start; A_end];

[T, T_dot] = ...
    compute_cheby(K,N,t);

[T_e, T_dot_e] = ...
    compute_cheby(K_e,N,t_e);

Phi = zeros(n,n*(N+1),K+1);
Phi_dot = zeros(n,n*(N+1),K+1);
for k = 1:K+1
    Phi(:,:,k) = kron(eye(n),T(:,k)');
    Phi_dot(:,:,k) = 2*kron(eye(n),T_dot(:,k)');
end

% Cost function
geo_cost_fnc =  @(vars) Geodesic_cost_tom(vars,w,n,...
    K,W,Phi,Phi_dot);

%Gradient function
geo_grad_fnc = @(vars) Geodesic_grad(vars,w,...
    K,N,n,T,W,dW,Phi,Phi_dot);

Name = 'Geodesic Problem';
geo_Prob = conAssign(geo_cost_fnc,geo_grad_fnc,[],[],...
                  [],[],Name,zeros(n*(N+1),1),[],0,...
                  Aeq,zeros(2*n,1),zeros(2*n,1),[],[],[],[],[],[]);
% geo_Prob.LineParam.MaxIter = 30;
% geo_Prob.LineParam.rho = 0.05;


end