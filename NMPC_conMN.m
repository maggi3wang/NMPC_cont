function c = NMPC_conMN(xu,D,n,m,N,f,B_full,Tp,P)
%Dynamics, and terminal

c = zeros(n*(N+1)+1,1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Tp)*D*xu(1:n*(N+1)) -...
    (NMPC_dyn(f,xu(1:n*(N+1)),n,N) + B_full*xu(n*(N+1)+1:end));

%% Terminal constraint
c(end) = xu(n*N+1:n*(N+1))'*P*xu(n*N+1:n*(N+1));

%% Initial RPI constraint

% x_bar = x_bar;
% c(end) = (xu(1:n)-x_bar)'*(xu(1:n)-x_bar);

end

