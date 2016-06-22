function conJ = NMPC_conJ(xu,Prob)%,D,n,m,N,df,B_full,Tp,P)
%Dynamics, and terminal
n = Prob.user.n;
N = Prob.user.N;
m = Prob.user.m;

global US_A;

conJ = zeros(n*(N+1)+2,(n+m)*(N+1));

%% Dynamics constraints

conJ(1:n*(N+1),1:n*(N+1)) = (2/Prob.user.Tp)*Prob.user.D - ...
            df_all(Prob.user.df,xu(1:n*(N+1)),n,N);

conJ(1:n*(N+1),n*(N+1)+1:end) = -Prob.user.B_full;

%% Initial RPI constraint

% conJ(end-1,1:n) = NaN;

M = (Prob.user.W(xu(1:n)))\eye(n);
conJ(end-1,1:n) = -2*US_A'*M;

% conJ(end-1,1:n) = -2*(Prob.user.M*(Prob.user.x_bar - xu(1:n)))';

%% Terminal constraint
conJ(end,n*N+1:n*(N+1)) = (2*Prob.user.P*(xu(n*N+1:n*(N+1))-Prob.user.x_eq))';



end

