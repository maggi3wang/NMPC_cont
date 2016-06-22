function J = brach_cond(xt,D,N,n)

X = xt(1:N+1);
Y = xt(N+2:2*(N+1));
theta = xt(2*(N+1)+1:end-1);
tau_f = xt(end);

J = zeros(2*(N+1),3*(N+1)+1);

J(1:N+1,1:N+1) = (2/tau_f)*D;
J(1:N+1,N+2:2*(N+1)) = (1/2)*sqrt(2)*diag(cos(theta)./sqrt(Y));
J(1:N+1,2*(N+1)+1:3*(N+1)) = sqrt(2)*diag(sqrt(Y).*sin(theta));
J(1:N+1,end) = -(2/(tau_f^2))*D*X;

J(N+2:2*(N+1),N+2:2*(N+1)) = (2/tau_f)*D - ...
                    (1/2)*sqrt(2)*diag(sin(theta)./sqrt(Y));
J(N+2:2*(N+1),2*(N+1)+1:3*(N+1)) = -sqrt(2)*diag(sqrt(Y).*cos(theta));
J(N+2:2*(N+1),end) = -(2/(tau_f^2))*D*Y;

                                           

end