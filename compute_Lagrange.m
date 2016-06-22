function L_e = compute_Lagrange(K,N,t,t_nodes)

L_e = ones(N+1,K+1); %up to order N

% %First get Chebyshev evals
% [T, T_dot] = compute_cheby(K,N,t);
% 
% %Only need TN_dot
% TN_dot = T_dot(N+1,:);

for k=0:N
%     if (j==0) || (j==N)
%         cj = 2;
%     else
%         cj = 1;
%     end
%     L_e(j+1,:) =  (((-1)^(j+1))/((N^2)*cj))*...
%                  (((1-t.^2)).*TN_dot./(t-t_nodes(j+1)));
    for m = 0:N
        if m~=k
            L_e(k+1,:) = L_e(k+1,:).*((t-t_nodes(m+1))/(t_nodes(k+1)-t_nodes(m+1)));
        end
    end
    
end


end
