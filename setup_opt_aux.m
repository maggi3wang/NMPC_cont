function aux_Prob = setup_opt_aux(K_e,Q_bar,m)

Name = 'Aux Problem';
aux_Prob = qpAssign(Q_bar,zeros(1,(K_e+1)*m),ones(K_e+1,(K_e+1)*m),...
         -Inf*ones((K_e+1)*m,1),zeros((K_e+1)*m,1),[],[],...
         zeros((K_e+1)*m,1),Name);

% aux_Prob = qpAssign(eye(m),zeros(1,m),ones(1,m),...
%          -Inf,0,[],[],...
%          0,Name);
       


end