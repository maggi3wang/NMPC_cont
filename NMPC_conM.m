function c = NMPC_conM(xu,D,n,N,f,B_full,Tf)
%Dynamics, and terminal

c = zeros(n*(N+1),1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Tf)*D*xu(1:n*(N+1)) -...
          (NMPC_dyn(f,xu(1:n*(N+1)),n,N) + B_full*xu(n*(N+1)+1:end));

end
