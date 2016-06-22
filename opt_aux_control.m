function u_aux = opt_aux_control(K_e,w_e,...
        edWfe,eWdfe,eB,eWe,X_dot,Q,Q_bar,lambda,m,n)

% X_dot_norms = norms(X_dot)';
% Q_bar = kron(diag(X_dot_norms),Q);

%%
% cvx_begin quiet
% variable K_vec(m*n,geo_prob_info.K + 1)
% variable K_vec(m*n*(geo_prob_info.K+1))
% J = 0;
% for k = 1:geo_prob_info.K+1
%     J = J + (1/2)*geo_prob_info.w(k)*...
%         norm(K_k(:,:,k),'fro')*X_dot_norms(k);
% end
% variable
% minimize (quad_form(K_vec,Q_bar)) 
% subject to
% for k = 1:geo_prob_info.K + 1
%     K_k = reshape(K_vec(:,k)',m,n);
%     -edWfe(k) + 2*eWdfe(k) + 2*(eB{k}*K_vec(:,k))'*We{k} <= -2*lambda*eWe(k);
%     -edWfe(k) + 2*eWdfe(k) + ...
%      2*(eB{k}*K_vec(1+(k-1)*m*n:k*m*n))'*We{k} <= ...
%     -2*lambda*eWe(k);
% end
% cvx_end

%%
% X_dot_stacked = reshape(X_dot,n*(geo_prob_info.K_e+1),1);
cvx_begin quiet
variable u_vec(m*(K_e+1)) 
% variable K(m,n,K_e+1)
minimize ( quad_form(u_vec,Q_bar) ) 
subject to
% for k = 1:K_e + 1
%     u_vec(:,k) == K(:,:,k)*X_dot(:,k);
% end
-edWfe + 2*eWdfe + 2*eB*u_vec <=-2*lambda*eWe;

cvx_end
%%          

u_aux = (1/2)*Q*u_vec;
% for k = 1:K_e+1
% %     K_k = reshape(K_vec(1+(k-1)*m*n:k*m*n)',m,n);
%     u_aux = u_aux + (1/2)*w_e(k)*u_vec(:,k);
% end

end