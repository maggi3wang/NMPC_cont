function GradObj = Geodesic_grad(vars,w,K,N,n,T,W_fnc,dW_fnc,Phi,Phi_dot)

GradObj = zeros(n*(N+1),1);
I = eye(n);

for k = 1:K+1 %0 ---> K
    x_k = Phi(:,:,k)*vars;
    x_dot_k = Phi_dot(:,:,k)*vars;
    
    W = W_fnc(x_k);
    M = W\eye(n);
    
    M_xdot = M*x_dot_k;
        
    W_dx = dW_fnc(x_k);
    
    GradObj = GradObj + w(k)*Phi_dot(:,:,k)'*M_xdot;
    
    for i = 1:n
        GradObj = GradObj+...
           -(1/2)*w(k)*(M_xdot'*W_dx{i}*M_xdot)*kron(I(:,i),T(:,k));
    end
end


return