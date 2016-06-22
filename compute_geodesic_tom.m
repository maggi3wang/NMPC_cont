function [X, X_dot,J_opt,exitflag] = compute_geodesic_tom(Prob,...
                                    n,N,start_p,end_p,T_e,T_dot_e,Aeq)


beq = [start_p; end_p];

vars_0 = zeros(n*(N+1),1);
for i = 1:n
    vars_0((i-1)*(N+1)+1:(i-1)*(N+1)+2) = [(start_p(i)+end_p(i))/2;
        -(start_p(i)-end_p(i))/2];
end

Prob = replace_A(Prob,Aeq,beq,beq); 
Prob = modify_x_0(Prob,vars_0);
Prob = ProbCheck(Prob, 'npsol');

Result = npsolTL(Prob);
                  
C_opt = (reshape(Result.x_k,N+1,n))';
X = C_opt*T_e;
X_dot = 2*C_opt*T_dot_e;

J_opt = Result.f_k;
exitflag = Result.Inform;

end



