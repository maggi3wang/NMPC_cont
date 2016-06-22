function [J,GradObj] = Geodesic_cost(vars,w,K,N,n,T,W_fnc,dW_fnc,Phi,Phi_dot)

GradObj = zeros(n*(N+1),1);

J = 0;

for k = 1:K+1 %0 ---> K
    x_k = Phi(:,:,k)*vars;
    x_dot_k = Phi_dot(:,:,k)*vars;
    
%     W = W0 + W1*x_k(1) + W2*(x_k(1)^2);
    W = W_fnc(x_k);
    M = W\eye(3);
    
    M_xdot = M*x_dot_k;
    J = J + (1/2)*w(k)*(x_dot_k'*M_xdot);
    
%     W_dx = W1 + 2*W2*x_k(1);
    W_dx = dW_fnc(x_k);
    
    GradObj = GradObj + w(k)*Phi_dot(:,:,k)'*M_xdot;
    
    GradObj = GradObj+...
            -(1/2)*w(k)*(M_xdot'*W_dx{1}*M_xdot)*...
                [T(:,k);zeros((n-1)*(N+1),1)];    
end



return