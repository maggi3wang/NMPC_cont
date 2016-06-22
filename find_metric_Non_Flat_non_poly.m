clear all; close all; clc;
yalmip('clear');

%% Allgower
box_size = 2;
ccm_eps = 0.1;
w_lower_eps = 0.1;
eps = 1e-3;

%lambda_range = 0.5:0.25:2;
% cond =1.000207519531250;
 lambda = 2 ;
 cond = 21.3984375;
bounds = nan(length(lambda_range),1);

% for ll = 1:length(lambda_range)
%     lambda = lambda_range(ll);
%     fprintf('lambda: %f ', lambda);
%     cond_l = 1;
%     cond_u = 50;
%     cond_best = cond_u;
%     found_feas = 0;
%     
%     while (cond_u-cond_l>eps)
%         cond = (cond_l+cond_u)/2;
        
        fprintf(' cond: %f ', cond);
%         lambda = 3;
        yalmip('clear');
        w_upper = sdpvar(1); 
        w_lower = sdpvar(1);
        
        
        x1 = sdpvar(1);
        x2 = sdpvar(1);
        
        x = [x1 x2]';
        delta = sdpvar(2,1);
        
        x_norm_sq = (x1^2) + (x2^2);
        sinh_term = x_norm_sq + (1/6)*(x_norm_sq^3);
        cosh_term = 1 + (1/2)*(x_norm_sq^2);
        
        f = [x2;
             -x1+x2*sinh_term];
        B = [0;
             1];
%         A = jacobian(f,[x1,x2]);
        A = [0, 1;
            -1+2*x1*x2*cosh_term,sinh_term+2*(x2^2)*cosh_term];
        
        Q_LQR = 10*eye(2); R_LQR = 2;
        [Klqr, P] = lqr(replace(A,x,zeros(2,1)),replace(B,x,zeros(2,1)),Q_LQR,R_LQR);
        Wlqr = inv(P);
                
        sol_order = 4;
        
        [W_11,c_11]= polynomial(x1,sol_order,0);
        [W_12,c_12]= polynomial(x1,sol_order,0);
        [W_22,c_22]= polynomial(x1,sol_order,0);
        
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
        s_order = 8;
        
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
                
        %W uniform bounds
        W_bounds = [w_upper>=w_lower,w_lower>=w_lower_eps];
        
        %Condition bound
        W_cond = [w_upper-cond*w_lower <= 0];
        
        %W pos def
        p_W_up = (w_upper*(delta'*delta)-delta'*W*delta)-[su_1,su_2,su_3,su_4]*g;
        
        p_W_low = (delta'*W*delta- w_lower*(delta'*delta))-[sl_1,sl_2,sl_3,sl_4]*g;
        
        %CCM condition
        R_CCM = -(-dW_f + A*W + W*A' - rho*(B*B') + 2*lambda*W);
        p_CCM = (delta'*R_CCM*delta - ccm_eps*(delta'*delta))-[s_ccm_1,s_ccm_2,s_ccm_3,s_ccm_4]*g;  
        
        %Assemble
        coeff_List = [c_rho;c_11;c_12;c_22;
                        c_us_1;c_us_2;c_us_3;c_us_4;c_ls_1;c_ls_2;c_ls_3;c_ls_4;...
                        c_ccm_1;c_ccm_2;c_ccm_3;c_ccm_4];
        
        constr_List = [ sos(p_W_low),sos(p_W_up),...
            sos(p_CCM),...
            sos(su_1), sos(su_2),sos(su_3), sos(su_4),...
            sos(sl_1), sos(sl_2),sos(sl_3), sos(sl_4),...
            W_bounds,W_cond,...
            sos(s_ccm_1),sos(s_ccm_2),sos(s_ccm_3),sos(s_ccm_4)];
        %         LQR_constr_W;LQR_constr_rho];
        
        options = sdpsettings('solver','mosek','verbose',1);
        SOS_soln = solvesos(constr_List,w_upper-w_lower,options,coeff_List);
        if (SOS_soln.problem == 0)
%             fprintf('feasible\n');
%             cond_u = cond;
%             cond_best = cond;
%             found_feas = 1;
            coeff_sol = clean(double(coeff_List),1e-3);
            W_sol = replace(W,coeff_List,coeff_sol);
            dW_x1 = jacobian(W(:),x1); dW_x1 = reshape(dW_x1,2,2);
            dW_x1_sol = replace(dW_x1,coeff_List,coeff_sol);
            rho_sol = replace(rho,coeff_List,coeff_sol);
        else
            fprintf('infeasible\n');
            cond_l = cond;
        end
%     end
%     fprintf('\n');
%     if (found_feas)
%         bounds(ll) = (sqrt(cond))/lambda;
%     end
% end
pause;
clc;

% figure()
% plot(lambda_range, bounds,'ro','markerfacecolor','g','markersize',20);
% grid on
% xlabel('\lambda');
% ylabel('$\|x^{*}-x\|/\bar{w}$','Interpreter','Latex');
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
% Display

disp('W')
sdisplay(W_sol)

disp('rho')
sdisplay(rho_sol)

disp('Bound achieved');
disp((1/lambda)*sqrt(value(cond)));

%% Convert into functions

W_cell = cell(2,2);
dW_x1_cell = cell(2,2);
for i = 1:2
    for j = 1:2
        W_ij = sdisplay(W_sol(i,j)); 
        dW_ij = sdisplay(dW_x1_sol(i,j));
        
        W_ij = strcat('@(x1,x2)',W_ij{1});
        dW_ij = strcat('@(x1)',dW_ij{1});
        
        W_cell{i,j} = str2func(W_ij);
        dW_x1_cell{i,j} = str2func(dW_ij);
    end
end
W_mat = @(x1,x2) [W_cell{1,1}(x1,x2), W_cell{1,2}(x1,x2);
                  W_cell{2,1}(x1,x2), W_cell{2,2}(x1,x2)];
dW_x1_mat = @(x1) [dW_x1_cell{1,1}(x1), dW_x1_cell{1,2}(x1);
                   dW_x1_cell{2,1}(x1), dW_x1_cell{2,2}(x1)];

rho = sdisplay(rho_sol);
rho = strcat('@(x1,x2)',rho{1});
rho = str2func(rho);

%% Plot

box_disp = 1.3;

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

f = @(x1,x2)[x2;
             -x1+x2*sinh(x1^2 + x2^2)];
B = [0;
     1];
A = @(x1,x2) [0, 1;
         -1+2*x1*x2*cosh(x1^2 + x2^2),sinh(x1^2 + x2^2)+2*(x2^2)*cosh(x1^2 + x2^2)];

% x_norm_sq = @(x1,x2) (x1^2) + (x2^2);
% sinh_term = @(x1,x2)  x_norm_sq(x1,x2) + (1/6)*(x_norm_sq(x1,x2)^3);
% cosh_term = @(x1,x2)  1 + (1/2)*(x_norm_sq(x1,x2)^2);
% 
% f = @(x1,x2)[x2;
%             -x1+x2*sinh_term(x1,x2)];
% B = [0;
%      1];
% A = @(x1,x2)[0, 1;
%            -1+2*x1*x2*cosh_term(x1,x2),sinh_term(x1,x2)+2*(x2^2)*cosh_term(x1,x2)];
 
R_val = zeros(size(X));

% (-dW_f + W*A' + A*W - rho*(B*B') + 2*lambda*W)
for i = 1:length(x2_range)    
    for j = 1:length(x1_range)
        x1 = x1_range(j); x2 = x2_range(i);
        
        L = (-dW_x1_mat(x1)*([1,0]*f(x1,x2)) + ...
                W_mat(x1,x2)*A(x1,x2)' + A(x1,x2)*W_mat(x1,x2) -...
                rho(x1,x2)*(B*B') + 2*lambda*W_mat(x1,x2));
            
        R_val(i,j) = max(eig(0.5*(L+L')));
    end
end

figure()
surface(X,Y,R_val);
grid on
hold on