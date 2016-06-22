function [X, X_dot,J_opt,exitflag] = compute_geodesic(prob_info,cost_fnc,start_p,end_p)


beq = [start_p; end_p];
n = prob_info.n;
N = prob_info.N;

vars_0 = zeros(n*(N+1),1);
for i = 1:n
    vars_0((i-1)*(N+1)+1:(i-1)*(N+1)+2) = [(start_p(i)+end_p(i))/2;
        -(start_p(i)-end_p(i))/2];
end

[vars_opt, J_opt, exitflag, solveroutput] = fmincon(cost_fnc,vars_0,[],[],...
    prob_info.Aeq,beq,[],[],[],prob_info.solver_options);

% J_opt = (norm(start_p-end_p))^2;
% exitflag = 1;
disp(vars_opt);
C_opt = (reshape(vars_opt,N+1,n))';
X = C_opt*prob_info.T_e;
X_dot = 2*C_opt*prob_info.T_dot_e;

end



