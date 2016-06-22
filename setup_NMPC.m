function [NMPC_Prob,Aeq,L_e,L_e_full] = ...
      setup_NMPC(n,m,f,B,df,x_L,u_L,N,Tp,delta,dt,P,alpha,W,geo_MPC,RPI_bound,x_eq)

%% Constants

%State bounds
% x_L = [-pi/4;-0.5;-0.4;-0.2];
x_U = -x_L;

%Control bounds
% u_L = -0.8*ones(m,1);
u_U = -u_L;

%State cost weighting
Q = diag([1,1, 15 ,1]);

%Control cost weighting
R = 10*eye(m);

%Terminal cost matrix (P,alpha)

%Order of approximation
% N = 10;

%Number of collocation points
K = N;

%MPC Solve horizon
% Tp = 10; %tau_p: [0, Tp]

%Interval
% delta = 1; 

%CGL nodes
[s_t,w] = clencurt(K); %t_t: [-1, 1] : <-> : [0, Tp]
s = fliplr(s_t); %t: [1, -1]

%% Final solution interpolation matrix

tau = 0:dt:delta;
s_e = (2*tau - Tp)/Tp; %[-1, s_delta]

tau_full = 0:dt:Tp; 
s_e_full = (2*tau_full - Tp)/Tp; %[-1, 1]

%Lagrange polynomial evaluation at the interpolation points
L_e = compute_Lagrange(length(s_e)-1,N,s_e,s_t);
L_e_full = compute_Lagrange(length(s_e_full)-1,N,s_e_full,s_t);
L_e_test = compute_Lagrange(length(s_t)-1,N,s_t,s_t);

%% Get Differentiation matrix

D = ChebyshevDiffMatrix(N,s); %arranged for forward time
D = kron(D,eye(n));

%% Variables

%State node values: [x_t0,...,x_tN]
%Control node values: [u_t0,...,u_tN]

% n_vars = (N+1)*(n+m);

%% Define problem

% x_eq = [0;0;0.2;0];
u_eq = zeros(m,1);

x_eq_all = kron(ones(N+1,1),x_eq);
u_eq_all = kron(ones(N+1,1),u_eq);

xu_eq = [x_eq_all;u_eq_all];

Q_bar = kron(diag([w(1:end-1);0]),Q); R_bar = kron(diag(w),R);
Q_tilde = Q_bar + kron(diag([zeros(N,1);(2/Tp)]),P);
% Q_tilde = Q_bar;
F = blkdiag(Q_tilde,R_bar);
F_pattern = sparse(F~=0);
    
B_full = kron(eye(N+1),B);

xu_L = [kron(ones(N+1,1),x_L);
        kron(ones(N+1,1),u_L)];
xu_U = [kron(ones(N+1,1),x_U);
        kron(ones(N+1,1),u_U)]; 
    
Aeq = [eye(n), zeros(n,N*n+(N+1)*m)];
b_L = zeros(n,1);
b_U = zeros(n,1);
    
MPC_cost = @(xu) (Tp/2)*(xu-xu_eq)'*F*(xu-xu_eq) ;
MPC_grad = @(xu) Tp*F*(xu-xu_eq);
% MPC_hess = Tp*F;
% nl_constr = @(xu) NMPC_con(xu,D,n,m,N,f,B_full,Tp,P);
% nl_conJ = @(xu) NMPC_conJ(xu,D,n,m,N,df,B_full,Tp,P);

c_L = zeros(n*(N+1)+2,1);
c_U = [zeros(n*(N+1),1);RPI_bound;alpha];

xu0 = zeros((n+m)*(N+1),1);

Name = 'NMPC';
NMPC_Prob = conAssign(MPC_cost,MPC_grad,[],[],...
            xu_L,xu_U,Name, xu0,...
            [], 0, [],[],[],...
            'NMPC_con','NMPC_conJ',[],[],...
            c_L,c_U,...
            [],[],[],[]);
        
% NMPC_Prob.CheckNaN = 1;
        
NMPC_Prob.user.x_bar = zeros(n,1);
NMPC_Prob.user.D = D;
NMPC_Prob.user.n = n;
NMPC_Prob.user.m = m;
NMPC_Prob.user.N = N;
NMPC_Prob.user.f = f;
NMPC_Prob.user.df = df;
NMPC_Prob.user.B_full = B_full;
NMPC_Prob.user.Tp = Tp;
NMPC_Prob.user.P = P;
% NMPC_Prob.user.M = M;
NMPC_Prob.user.W = W;
NMPC_Prob.user.geo_MPC = geo_MPC;
NMPC_Prob.user.x_eq = x_eq;
%% Toy Brach Problem

% Aeq = [eye(n), zeros(n,n*N+N+1+1);
%        zeros(1,n*N),1,0,zeros(1,N+1+1)];
%    
% b_L = [zeros(2,1);0.5];
% b_U = b_L;  
% 
% nl_constr = @(xt) brach_con(xt,D,N,n);
% % nl_conJ = @(xt) brach_cond(xt,D,N,n);
% n_vars = 3*(N+1)+1;
% x0 = [kron([zeros(N,1);0.5],[1;0])+...
%       kron((linspace(0,1,N+1)+0.1)',[0;1]);
%       zeros(N+1,1);0.4];
% x_L = [zeros(2*(N+1),1);-pi*ones(N+1,1);0];
% x_U = [5*ones(2*(N+1),1);pi*ones(N+1,1);2];
% 
% Name = 'Brach';
% NMPC_Prob = conAssign(@(xt) xt(end),[],[],[],...
%             x_L,x_U,Name, x0,...
%             [], 0, Aeq, b_L,b_U,...
%             nl_constr, [],[],[],...
%             zeros(2*(N+1),1),zeros(2*(N+1),1),...
%             [],[],[],[]);






end