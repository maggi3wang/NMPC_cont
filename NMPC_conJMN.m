function conJ = NMPC_conJMN(xu,D,n,m,N,df,B_full,Tp,P)
%Dynamics, and terminal

conJ = zeros(n*(N+1)+1,(n+m)*(N+1));

%% Dynamics constraints

conJ(1:n*(N+1),1:n*(N+1)) = (2/Tp)*D - ...
            df_all(df,xu(1:n*(N+1)),n,N);

conJ(1:n*(N+1),n*(N+1)+1:end) = -B_full;

%% Terminal constraint
conJ(end,n*N+1:n*(N+1)) = (2*P*xu(n*N+1:n*(N+1)))';

end

