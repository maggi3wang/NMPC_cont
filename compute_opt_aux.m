function [u_aux,e_rate,converged] = compute_opt_aux(Prob,Q,w_e,K_e,E,X,X_dot,W,dW,f,df,B,u_nom,lambda,n,m)


edWfe = zeros(K_e+1,1);
eWdfe = zeros(K_e+1,1);
eWe = E*ones(K_e+1,1);
eB = zeros(K_e+1,(K_e+1)*m);

for k = 1:K_e+1
   
    W_geo = W(X(:,k));
    Eta = W_geo\X_dot(:,k);
    
    dWf = zeros(n,n);
    dW_geo = dW(X(:,k));
    f_geo = f(X(:,k));
    for i = 1:n
        dWf = dWf + dW_geo{i}*f_geo(i);
    end
    edWfe(k) = Eta'*dWf*Eta;
    We = W_geo*Eta;
    eWdfe(k) = (We')*df(X(:,k))'*Eta;
    eB(k,1+(k-1)*m:k*m) = Eta'*B;
%     eWe(k) = Eta'*We;
end


Prob = replace_A(Prob,2*eB,-Inf*ones((K_e+1)*m,1),-2*lambda*eWe + edWfe - 2*eWdfe);

Prob = ProbCheck(Prob, 'qpopt');
Result = qpoptTL(Prob);

delta_u = Result.x_k;

for k = 1:K_e+1
    W_geo = W(X(:,k));
    Eta = W_geo\X_dot(:,k);
    if norm(Eta'*B)<1e-3
%         disp('found');
        delta_u(1+(k-1)*m:k*m) = zeros(m,1);
    end
end

e_rate = (1/2)*w_e'*(-edWfe + 2*eWdfe + 2*eB*delta_u);
% if e_rate >0
% %     keyboard;
%     u_aux = 0;
% else
    u_aux = (1/2)*Q*delta_u;
% end

if norm(u_aux)>1
%     keyboard;
end

%%

% k = K_e + 1;
% A = 2*X_dot(:,k)'*(W(X(:,k))\B);
% l_b = -Inf;
% u_b = -2*lambda*E + 2*X_dot(:,1)'*(W(X(:,1))\(f(X(:,1)) + B*u_nom)) -...
%              2*X_dot(:,k)'*(W(X(:,k))\(f(X(:,k)) + B*u_nom));
% Prob = replace_A(Prob,A,l_b,u_b);

% e_rate =     2*X_dot(:,k)'*(W(X(:,k))\(f(X(:,k)) + B*Result.x_k)) -...
%              2*X_dot(:,1)'*(W(X(:,1))\(f(X(:,1)) + B*u_nom));

% u_aux = Result.x_k;

%%
converged = Result.Inform;

end



