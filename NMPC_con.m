function c = NMPC_con(xu,Prob)%,D,n,m,N,f,B_full,Tp,P)
%Dynamics, and terminal

n = Prob.user.n;
N = Prob.user.N;

geo = Prob.user.geo_MPC;

global US_A;

c = zeros(n*(N+1)+2,1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Prob.user.Tp)*Prob.user.D*xu(1:n*(N+1)) -...
    (NMPC_dyn(Prob.user.f,xu(1:n*(N+1)),n,N) + Prob.user.B_full*xu(n*(N+1)+1:end));

%% Initial RPI constraint

[~, X_dot,J_opt,~] = compute_geodesic_tom(geo.geo_Prob,...
    n,geo.geodesic_N,xu(1:n),Prob.user.x_bar,geo.T_e,geo.T_dot_e,geo.geo_Aeq);
US_A = X_dot(:,1);

c(end-1) = J_opt;

% c(end-1) = (Prob.user.x_bar- xu(1:n))'*(Prob.user.M)*(Prob.user.x_bar-xu(1:n));

%% Terminal constraint
c(end) = (xu(n*N+1:n*(N+1))-Prob.user.x_eq)'*Prob.user.P*(xu(n*N+1:n*(N+1))-Prob.user.x_eq);



end

