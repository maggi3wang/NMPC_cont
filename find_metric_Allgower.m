clear all; close all; clc;
yalmip('clear');

%% Allgower

lambda_range = 1.742857142857143;
% cond =1.000207519531250;
% lambda_range = linspace(0.8,2,15);
bounds = nan(15,1);
eps = 1e-3;
for ll = 1:length(lambda_range)
    lambda = lambda_range(ll);
    fprintf('lambda: %f ', lambda);
    cond_l = 1;
    cond_u = 5;
    cond_best = cond_u;
    
    while (cond_u-cond_l>eps)
        cond = (cond_l+cond_u)/2;
        
        fprintf(' cond: %f ', cond);
        
        yalmip('clear');
        w_upper = sdpvar(1); w_lower = sdpvar(1);
        
        x1 = sdpvar(1);
        x2 = sdpvar(1);
        
        x = [x1 x2]';
        
        f = [-1*x1 + 2*x2;
            -3*x1 + 4*x2 - 0.25*(x2^3)];
        B = [0.5;-2];
        A = jacobian(f,x);
        
        A0 = replace(A,x,zeros(2,1));
        Q_LQR = 10*eye(2); R_LQR = 2;
        [Klqr, P] = lqr(A0,B,Q_LQR,R_LQR);
        Wlqr = inv(P);
        
        
        delta = sdpvar(2,1);
        
        [W_11,c_11]= polynomial(x,4,0);
        [W_12,c_12]= polynomial(x,4,0);
        [W_22,c_22]= polynomial(x,4,0);
        
        [rho,c_rho] = polynomial(x,4,0);
        
        W = [W_11,W_12;
            W_12,W_22];
        
        dW_f = jacobian(W(:),x)*f;
        dW_f = reshape(dW_f,2,2);
        
        % Enforce LQR at origin
        LQR_constr_W   = [replace(W,x,0) == Wlqr];
        LQR_constr_rho = [replace(rho,x,0) == 2*inv(R_LQR)];
        
       
        %Box constraints
        g = [x1+5; 5-x1; x2+5; 5-x2];
        
        %Box Lagrange functions
        
        %Upper and Lower definiteness
        [su_1,c_us_1] = polynomial([x;delta],4,0);
        [su_2,c_us_2] = polynomial([x;delta],4,0);
        [su_3,c_us_3] = polynomial([x;delta],4,0);
        [su_4,c_us_4] = polynomial([x;delta],4,0);
        [sl_1,c_ls_1] = polynomial([x;delta],4,0);
        [sl_2,c_ls_2] = polynomial([x;delta],4,0);
        [sl_3,c_ls_3] = polynomial([x;delta],4,0);
        [sl_4,c_ls_4] = polynomial([x;delta],4,0);
        
        %Killing Field 
        [seq_11_1,ceq_11_1] = polynomial(x,4,0);
        [seq_11_2,ceq_11_2] = polynomial(x,4,0);
        [seq_11_3,ceq_11_3] = polynomial(x,4,0);
        [seq_11_4,ceq_11_4] = polynomial(x,4,0);
        [seq_11_5,ceq_11_5] = polynomial(x,4,0);
        [seq_11_6,ceq_11_6] = polynomial(x,4,0);
        [seq_11_7,ceq_11_7] = polynomial(x,4,0);
        [seq_11_8,ceq_11_8] = polynomial(x,4,0);
        
        [seq_12_1,ceq_12_1] = polynomial(x,4,0);
        [seq_12_2,ceq_12_2] = polynomial(x,4,0);
        [seq_12_3,ceq_12_3] = polynomial(x,4,0);
        [seq_12_4,ceq_12_4] = polynomial(x,4,0);
        [seq_12_5,ceq_12_5] = polynomial(x,4,0);
        [seq_12_6,ceq_12_6] = polynomial(x,4,0);
        [seq_12_7,ceq_12_7] = polynomial(x,4,0);
        [seq_12_8,ceq_12_8] = polynomial(x,4,0);
        
        [seq_22_1,ceq_22_1] = polynomial(x,4,0);
        [seq_22_2,ceq_22_2] = polynomial(x,4,0);
        [seq_22_3,ceq_22_3] = polynomial(x,4,0);
        [seq_22_4,ceq_22_4] = polynomial(x,4,0);
        [seq_22_5,ceq_22_5] = polynomial(x,4,0);
        [seq_22_6,ceq_22_6] = polynomial(x,4,0);
        [seq_22_7,ceq_22_7] = polynomial(x,4,0);
        [seq_22_8,ceq_22_8] = polynomial(x,4,0);
        
        %CCM condition
        [s_ccm_1,c_ccm_1] = polynomial([x;delta],4,0);
        [s_ccm_2,c_ccm_2] = polynomial([x;delta],4,0);
        [s_ccm_3,c_ccm_3] = polynomial([x;delta],4,0);
        [s_ccm_4,c_ccm_4] = polynomial([x;delta],4,0);
        
        %CCM condition
        R_CCM = -(-dW_f + A*W + W*A' - rho*(B*B') + 2*lambda*W);
        p_CCM = delta'*R_CCM*delta - [s_ccm_1,s_ccm_2,s_ccm_3,s_ccm_4]*g;        
        
        %Killing field condition
        dW_B = jacobian(W(:),x)*B;
        dW_B = reshape(dW_B,2,2);
        killing_11 = dW_B(1,1) - [seq_11_1,seq_11_2,seq_11_3,seq_11_4]*g;
        killing_11_2 = (-dW_B(1,1)) - [seq_11_5,seq_11_6,seq_11_7,seq_11_8]*g;
        
        killing_12 = dW_B(1,2) - [seq_12_1,seq_12_2,seq_12_3,seq_12_4]*g;
        killing_12_2 = (-dW_B(1,2)) - [seq_12_5,seq_12_6,seq_12_7,seq_12_8]*g;
        
        killing_22 = dW_B(2,2) - [seq_22_1,seq_22_2,seq_22_3,seq_22_4]*g;
        killing_22_2 = (-dW_B(2,2)) - [seq_22_5,seq_22_6,seq_22_7,seq_22_8]*g;
        
        %W uniform bounds
        W_bounds = [w_lower>=0.0035; w_upper>=w_lower];
        
        %Condition bound
        W_cond = [w_upper-cond*w_lower <= 0];
        
        %W pos def
        p_W_up = (w_upper*(delta'*delta)-delta'*W*delta) - ...
            [su_1,su_2,su_3,su_4]*g;
        
        p_W_low = (delta'*W*delta - w_lower*(delta'*delta))-...
            [sl_1,sl_2,sl_3,sl_4]*g;
        
        coeff_List = [c_rho;c_11;c_12;c_22;
            c_us_1;c_us_2;c_us_3;c_us_4;c_ls_1;c_ls_2;c_ls_3;c_ls_4;
            ceq_11_1;ceq_11_2;ceq_11_3;ceq_11_4;ceq_11_5;ceq_11_6;ceq_11_7;ceq_11_8;
            ceq_12_1;ceq_12_2;ceq_12_3;ceq_12_4;ceq_12_5;ceq_12_6;ceq_12_7;ceq_12_8;
            ceq_22_1;ceq_22_2;ceq_22_3;ceq_22_4;ceq_22_5;ceq_22_6;ceq_22_7;ceq_22_8;
            c_ccm_1;c_ccm_2;c_ccm_3;c_ccm_4];
        
        constr_List = [ sos(p_W_low);sos(p_W_up);
            sos(p_CCM);
            sos(su_1); sos(su_2);
            sos(su_3); sos(su_4);
            sos(sl_1); sos(sl_2);
            sos(sl_3); sos(sl_4);
            W_bounds;W_cond;
            sos(seq_11_1);sos(seq_11_2);sos(seq_11_3);sos(seq_11_4);
            sos(seq_11_5);sos(seq_11_6);sos(seq_11_7);sos(seq_11_8);
            sos(seq_12_1);sos(seq_12_2);sos(seq_12_3);sos(seq_12_4);
            sos(seq_12_5);sos(seq_12_6);sos(seq_12_7);sos(seq_12_8);
            sos(seq_22_1);sos(seq_22_2);sos(seq_22_3);sos(seq_22_4);
            sos(seq_22_5);sos(seq_22_6);sos(seq_22_7);sos(seq_22_8);
            sos(killing_11); sos(killing_11_2);
            sos(killing_12); sos(killing_12_2);
            sos(killing_22); sos(killing_22_2);
            sos(s_ccm_1); sos(s_ccm_2); sos(s_ccm_3); sos(s_ccm_4)];
        %         LQR_constr_W;LQR_constr_rho];
        
        options = sdpsettings('solver','mosek','verbose',0);
        SOS_soln = solvesos(constr_List,[],options,coeff_List);
        if (SOS_soln.problem == 0)
            fprintf('feasible\n');
            cond_u = cond;
            cond_best = cond;
            coeff_sol = clean(double(coeff_List),1e-3);
            W_sol = replace(W,coeff_List,coeff_sol);
            dW_f = jacobian(W(:),x)*f; dW_f = reshape(dW_f,2,2);
            dW_f_sol = replace(dW_f,coeff_List,coeff_sol);
            rho_sol = replace(rho,coeff_List,coeff_sol);
        else
            cond_l = cond;
        end
    end
    fprintf('\n');
    bounds(ll) = (sqrt(cond))/lambda;
end
pause;
clc;

figure()
plot(lambda_range, bounds,'ro','markerfacecolor','g','markersize',20);
grid on
xlabel('\lambda');
ylabel('$\|x^{*}-x\|/\bar{w}$','Interpreter','Latex');
% title('Robustness optimization');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
% Display

disp('W')
sdisplay(W_sol)

disp('rho')
sdisplay(rho_sol)

disp('Bound achieved');
disp((1/lambda)*sqrt(cond));

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

x1_range = -5:0.2:5;
x2_range = -5:0.2:5;
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

f = @(x1,x2) [-1*x1 + 2*x2;
             -3*x1 + 4*x2 - 0.25*(x2^3)];
B = [0.5;-2];
A = @(x1,x2) [-1, 2;
              -3, 4-0.75*(x2^2)];
R_val = zeros(size(X));

% (-dW_f + A*W + W*A' - rho*(B*B') + 2*lambda*W)
for i = 1:length(x2_range)    
    for j = 1:length(x1_range)
        x1 = x1_range(j); x2 = x2_range(i);
        
        L = (-dW_f_mat(x1,x2) + ...
                A(x1,x2)*W_mat(x1,x2) + W_mat(x1,x2)*A(x1,x2)' -...
                rho(x1,x2)*(B*B') + 2*lambda*W_mat(x1,2));
        R_val(i,j) = max(eig(L));
    end
end

figure()
surface(X,Y,R_val);
grid on
hold on