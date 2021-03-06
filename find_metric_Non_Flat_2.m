clear all; close all; clc;
yalmip('clear');

%% Allgower
box_size = 5;
% lambda_range = 1.742857142857143;
% cond =1.000207519531250;
% lambda_range = linspace(0.8,2,15);
% bounds = nan(15,1);
% eps = 1e-3;
% for ll = 1:length(lambda_range)
%     lambda = lambda_range(ll);
%     fprintf('lambda: %f ', lambda);
%     cond_l = 1;
%     cond_u = 5;
%     cond_best = cond_u;
%     
%     while (cond_u-cond_l>eps)
%         cond = (cond_l+cond_u)/2;
%         
%         fprintf(' cond: %f ', cond);
        lambda = 20;
        yalmip('clear');
        w_upper = sdpvar(1); 
        w_lower = sdpvar(1);
%         ccm_eps = 0.1;
        w_lower_eps = 0.05;
        
        x1 = sdpvar(1);
        x2 = sdpvar(1);
        
        x = [x1 x2]';
        delta = sdpvar(2,1);
        
        f = [-50*x1 - 10*(x1^2);
             50*x1 - 100*x2];
        B = [10-x1;
             -x2];
        A = jacobian(f,[x1,x2]);
        
        Q_LQR = 10*eye(2); R_LQR = 2;
        [Klqr, P] = lqr(replace(A,x,zeros(2,1)),replace(B,x,zeros(2,1)),Q_LQR,R_LQR);
        Wlqr = inv(P);
                
        sol_order = 4;
        
        [W_11,c_11]= polynomial(x,sol_order,0);
        [W_12,c_12]= polynomial(x,sol_order,0);
        [W_22,c_22]= polynomial(x,sol_order,0);
        
        [rho,c_rho] = polynomial(x,sol_order,0);
        
        W = [W_11,W_12;
             W_12,W_22];
        
        dW_f = jacobian(W(:),x')*f;
        dW_f = reshape(dW_f,2,2);
        
        % Enforce LQR at origin
        LQR_constr_W   = [replace(W,x,0) == Wlqr];
        LQR_constr_rho = [replace(rho,x,0) == 2*inv(R_LQR)];
        
       
        %Box constraints
        g = [x1+box_size; box_size-x1; x2+box_size; box_size-x2];
        
        %Box Lagrange functions
        s_order = 6;
        
        %Upper and Lower definiteness
        [su_1,c_us_1] = polynomial([x;delta],s_order,0);
        [su_2,c_us_2] = polynomial([x;delta],s_order,0);
        [su_3,c_us_3] = polynomial([x;delta],s_order,0);
        [su_4,c_us_4] = polynomial([x;delta],s_order,0);
        [sl_1,c_ls_1] = polynomial([x;delta],s_order,0);
        [sl_2,c_ls_2] = polynomial([x;delta],s_order,0);
        [sl_3,c_ls_3] = polynomial([x;delta],s_order,0);
        [sl_4,c_ls_4] = polynomial([x;delta],s_order,0);
                
        %CCM condition
        [s_ccm_1,c_ccm_1] = polynomial([x;delta],s_order,0);
        [s_ccm_2,c_ccm_2] = polynomial([x;delta],s_order,0);
        [s_ccm_3,c_ccm_3] = polynomial([x;delta],s_order,0);
        [s_ccm_4,c_ccm_4] = polynomial([x;delta],s_order,0);    
        
        %Killing field condition
        dW_B = jacobian(W(:),[x1,x2])*B;
        dW_B = reshape(dW_B,2,2);
        
        W_db = -W;
            
        killing_11_1 =   dW_B(1,1)-2*W_db(1,1);
        killing_11_2 = -(dW_B(1,1)-2*W_db(1,1));
        
        killing_12_1 =   dW_B(1,2)-W_db(1,2)-W_db(2,1);
        killing_12_2 = -(dW_B(1,2)-W_db(1,2)-W_db(2,1));
        
        killing_22_1 =   dW_B(2,2)-2*W_db(2,2);
        killing_22_2 = -(dW_B(2,2)-2*W_db(2,2));
        
        %W uniform bounds
        W_bounds = [w_upper>=w_lower,w_lower>=w_lower_eps];
        
        %Condition bound
%         W_cond = [w_upper-cond*w_lower <= 0];
        
        %W pos def
        p_W_up = (w_upper*(delta'*delta)-delta'*W*delta)-[su_1,su_2,su_3,su_4]*g;
        
        p_W_low = (delta'*W*delta- w_lower*(delta'*delta))-[sl_1,sl_2,sl_3,sl_4]*g;
        
        %CCM condition
        R_CCM = -(-dW_f + A*W + W*A' - rho*(B*B') + 2*lambda*W);
        p_CCM = delta'*R_CCM*delta-[s_ccm_1,s_ccm_2,s_ccm_3,s_ccm_4]*g;  
        
        %Assemble
        coeff_List = [c_rho;c_11;c_12;c_22;
                        c_us_1;c_us_2;c_us_3;c_us_4;c_ls_1;c_ls_2;c_ls_3;c_ls_4;...
                        c_ccm_1;c_ccm_2;c_ccm_3;c_ccm_4];
        
        constr_List = [ sos(p_W_low),sos(p_W_up),...
            sos(p_CCM),...
            sos(su_1), sos(su_2),sos(su_3), sos(su_4),...
            sos(sl_1), sos(sl_2),sos(sl_3), sos(sl_4),...
            W_bounds,...%W_cond;
            sos(killing_11_1), sos(killing_11_2),...
            sos(killing_12_1), sos(killing_12_2),...
            sos(killing_22_1), sos(killing_22_2),...
            sos(s_ccm_1),sos(s_ccm_2),sos(s_ccm_3),sos(s_ccm_4)];
        %         LQR_constr_W;LQR_constr_rho];
        
        options = sdpsettings('solver','mosek','verbose',1);
        SOS_soln = solvesos(constr_List,w_upper-w_lower,options,coeff_List);
%         if (SOS_soln.problem == 0)
%             fprintf('feasible\n');
%             cond_u = cond;
%             cond_best = cond;
            coeff_sol = clean(double(coeff_List),1e-3);
            W_sol = replace(W,coeff_List,coeff_sol);
            dW_f = jacobian(W(:),x)*f; dW_f = reshape(dW_f,2,2);
            dW_f_sol = replace(dW_f,coeff_List,coeff_sol);
            rho_sol = replace(rho,coeff_List,coeff_sol);
%         else
%             cond_l = cond;
%         end
%     end
%     fprintf('\n');
%     bounds(ll) = (sqrt(cond))/lambda;
% end
pause;
clc;

% figure()
% plot(lambda_range, bounds,'ro','markerfacecolor','g','markersize',20);
% grid on
% xlabel('\lambda');
% ylabel('$\|x^{*}-x\|/\bar{w}$','Interpreter','Latex');
% % title('Robustness optimization');
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
% Display

disp('W')
sdisplay(W_sol)

disp('rho')
sdisplay(rho_sol)

disp('Bound achieved');
disp((1/lambda)*sqrt(value(w_upper/w_lower)));

%% Convert into functions
W_cell = cell(2,2);
dW_f_cell = cell(2,2);
for i = 1:2
    for j = 1:2
        W_ij = sdisplay(W_sol(i,j)); 
        dW_ij = sdisplay(dW_f_sol(i,j));
        
        W_ij = strcat('@(x1,x2)',W_ij{1});
        dW_ij = strcat('@(x1,x2)',dW_ij{1});
        
        W_cell{i,j} = str2func(W_ij);
        dW_f_cell{i,j} = str2func(dW_ij);
    end
end
W_mat = @(x1,x2) [W_cell{1,1}(x1,x2), W_cell{1,2}(x1,x2);
                  W_cell{2,1}(x1,x2), W_cell{2,2}(x1,x2)];
dW_f_mat = @(x1,x2) [dW_f_cell{1,1}(x1,x2), dW_f_cell{1,2}(x1,x2);
                   dW_f_cell{2,1}(x1,x2), dW_f_cell{2,2}(x1,x2)];

rho = sdisplay(rho_sol);
rho = strcat('@(x1,x2)',rho{1});
rho = str2func(rho);

box_disp = box_size;

x1_range = -box_disp:0.2:box_disp;
x2_range = -box_disp:0.2:box_disp;
[X,Y] = meshgrid(x1_range,x2_range);
W_max = zeros(size(X));
W_min = zeros(size(X));

for i = 1:length(x2_range)
    for j = 1:length(x1_range)
        l = eig(W_mat(x1_range(j),x2_range(i)));
        W_max(i,j) = max(l);
        W_min(i,j) = min(l);
    end
end

figure()
surface(X,Y,W_min);
hold on
surface(X,Y,W_max);
grid on

% Check Contraction condition

f = @(x1,x2)[-50*x1 - 10*(x1^2);
             50*x1 - 100*x2];
B = @(x1,x2)[10-x1;
             -x2];
A = @(x1,x2)[-50-20*x1, 0;
              50, -100];
R_val = zeros(size(X));

% (-dW_f + A*W + W*A' - rho*(B*B') + 2*lambda*W)
for i = 1:length(x2_range)    
    for j = 1:length(x1_range)
        x1 = x1_range(j); x2 = x2_range(i);
        
        L = (-dW_f_mat(x1,x2) + ...
                W_mat(x1,x2)*A(x1,x2)' + A(x1,x2)*W_mat(x1,x2) -...
                rho(x1,x2)*(B(x1,x2)*B(x1,x2)') + 2*lambda*W_mat(x1,x2));
            
        R_val(i,j) = max(eig(0.5*(L+L')));
    end
end

figure()
surface(X,Y,R_val);
grid on
hold on