function [NMPC_state,NMPC_ctrl,x_nom,u_nom,converged] = ...
    compute_MNMPC(Prob,Aeq,start_p,n,m,N,L_e)

Prob = replace_A(Prob,Aeq,start_p,start_p);
In = eye(n);
x0 = zeros((N+1)*n,1);
for i = 1:n
    x0 = x0 + kron(linspace(start_p(i),0,N+1), In(i,:))';
end
Prob = modify_x_0(Prob,x0,1:(N+1)*n);

Prob = ProbCheck(Prob, 'snopt');
Result = snoptTL(Prob);

converged = ~Result.ExitFlag;

%Compute trajectories
NMPC_state = zeros(size(L_e,2),n);
x_nom = zeros(N+1,n);
for i = 1:n
    c = Result.x_k(i:n:n*(N+1)-(n-i))';
    NMPC_state(:,i) = (c*L_e)';
    x_nom(:,i) = c';
end

NMPC_ctrl = zeros(size(L_e,2),m);
u_nom = zeros(N+1,m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-(m-j))';
    NMPC_ctrl(:,j) = (c*L_e)';
    u_nom(:,j) = c';
end

end