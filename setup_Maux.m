function Maux_prob = setup_Maux(n,m,B,N)

%% Constants

%State bounds
x_L = -7*ones(n,1);
x_U = -x_L;

%Control bounds
u_L = -2*ones(m,1);
u_U = -u_L;

%State cost weighting
Q = 0.5*eye(n);

%Control cost weighting
R = 1*eye(m);

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

Maux_prob.s_t = s_t;
Maux_prob.w = w;

%% Final solution interpolation matrix

%% Get Differentiation matrix

D = ChebyshevDiffMatrix(N,s);
Maux_prob.D = kron(D,eye(n));

%% Variables

%State node values: [x_t0,...,x_tN]
%Control node values: [u_t0,...,u_tN]

% n_vars = (N+1)*(n+m);

%% Define problem

Maux_prob.Q_bar = kron(diag(w),Q); 
Maux_prob.R_bar = kron(diag(w),R);
    
Maux_prob.B_full = kron(eye(N+1),B);

Maux_prob.xu_L = [kron(ones(N+1,1),x_L);
                  kron(ones(N+1,1),u_L)];
Maux_prob.xu_U = [kron(ones(N+1,1),x_U);
                  kron(ones(N+1,1),u_U)]; 

% Maux_prob.xu_L = [-Inf*ones((N+1)*n,1);
%                    kron(ones(N+1,1),u_L)];
% Maux_prob.xu_U = [Inf*ones((N+1)*n,1);
%                    kron(ones(N+1,1),u_U)];    
Maux_prob.Aeq = [eye(n), zeros(n,N*n+(N+1)*m)];



end