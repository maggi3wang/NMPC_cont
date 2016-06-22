function [NMPC_state,NMPC_ctrl,converged] = ...
    compute_NMPC(Prob,Aeq,start_p,mpc_state,state_constr,x_eq,n,m,N,L_e)

% Prob = replace_A(Prob,Aeq,start_p,start_p);

for i = 1:n
    if abs(mpc_state(i)) > abs(state_constr(i));
        mpc_state(i) = 0.98*sign(mpc_state(i))*abs(state_constr(i));
    end
end

In = eye(n);
x0 = zeros((N+1)*n,1);
for i = 1:n
    x0 = x0 + kron(linspace(mpc_state(i),x_eq(i),N+1), In(i,:))';
end
Prob = modify_x_0(Prob,x0,1:(N+1)*n);
Prob.user.x_bar = start_p;

Prob = ProbCheck(Prob, 'snopt');
Result = snoptTL(Prob);

converged = Result.Inform;

%Compute trajectories
NMPC_state = zeros(size(L_e,2),n);
for i = 1:n
    c = Result.x_k(i:n:n*(N+1)-(n-i))';
    NMPC_state(:,i) = (c*L_e)';
end

NMPC_ctrl = zeros(size(L_e,2),m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-(m-j))';
    NMPC_ctrl(:,j) = (c*L_e)';
%     if (sum(abs(NMPC_ctrl(:,j))>1.9)~=0)
%         keyboard;
%     end
end

end