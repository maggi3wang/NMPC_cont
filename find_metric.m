clear all; close all; clc;
yalmip('clear');

%% Slotine Dynamics

% y_w_upper_0 = log(1.561470885729414e+02);
% y_w_lower_0 = log(0.009979);
% y_lambda_0 = log(0.731623);

% lambda_range = linspace(0.1,1,15);
lambda_range = 0.6;
cond_l = 0.1;
cond_u = 500;
eps = 1e-3;
cond_best = cond_u;
% 
while (cond_u-cond_l>eps)
    cond = (cond_l+cond_u)/2;
%         cond =1.393188476562500;
    fprintf('cond: %f \n', cond);
    lambda = lambda_range(1);

    yalmip('clear');
    w_upper = sdpvar(1); w_lower = sdpvar(1);
    
    x1 = sdpvar(1);
    x2 = sdpvar(1);
    x3 = sdpvar(1);
    
    x = [x1 x2 x3]';
    
    f = [-x1+x3;
        x1^2-x2-2*x1*x3+x3;
        -x2];
    B = [0;0;1];
    A = jacobian(f,x);
    
    A0 = replace(A,x,zeros(3,1));
    Q_LQR = 50*eye(3); R_LQR = 1;
    [Klqr, P] = lqr(A0,B,Q_LQR,R_LQR);
    Wlqr = inv(P);
    
    w_states = [x1; x2];
    delta = sdpvar(3,1);
    
    [W_11,c_11]= polynomial(w_states,4,0);
    [W_12,c_12]= polynomial(w_states,4,0);
    [W_13,c_13]= polynomial(w_states,4,0);
    [W_22,c_22]= polynomial(w_states,4,0);
    [W_23,c_23]= polynomial(w_states,4,0);
    [W_33,c_33]= polynomial(w_states,4,0);
    
    [rho,c_rho] = polynomial(w_states,4,0);
    
    W = [W_11,W_12,W_13;
        W_12,W_22,W_23;
        W_13,W_23,W_33];
    
    W_dot = jacobian(W(:),x)*f;
    W_dot = reshape(W_dot,3,3);
    
    % Enforce LQR at origin
    LQR_constr_W   = [replace(W,x,0) == Wlqr];
    LQR_constr_rho = [replace(rho,x,0) == 2*inv(R_LQR)];
    
    %CCM condition
    R_CCM = -(-W_dot + A*W + W*A' - rho*(B*B') + 2*lambda*W);
    p_CCM = delta'*R_CCM*delta;
    
    %Box constraints
    g = [x1+4; 4-x1; x2+4; 4-x2];
%     
    %Box Lagrange functions
    [su_1,c_us_1] = polynomial([w_states;delta],4,0);
    [su_2,c_us_2] = polynomial([w_states;delta],4,0);
    [su_3,c_us_3] = polynomial([w_states;delta],4,0);
    [su_4,c_us_4] = polynomial([w_states;delta],4,0);
    [sl_1,c_ls_1] = polynomial([w_states;delta],4,0);
    [sl_2,c_ls_2] = polynomial([w_states;delta],4,0);
    [sl_3,c_ls_3] = polynomial([w_states;delta],4,0);
    [sl_4,c_ls_4] = polynomial([w_states;delta],4,0);
    
    
    %W uniform bounds
    W_bounds = [w_lower>=0.02; w_upper>=0.0035];
%     
%     %Condition bound
%     W_cond = [w_upper-cond*w_lower <= 0];
%     
    %W pos def
    p_W_up = (w_upper*(delta'*delta)-delta'*W*delta) - ...
        [su_1,su_2,su_3,su_4]*g;
    
    p_W_low = (delta'*W*delta - w_lower*(delta'*delta))-...
        [sl_1,sl_2,sl_3,sl_4]*g;
    
    
    coeff_List = [c_rho;c_11;c_12;c_13;c_22;c_23;c_33;
        c_us_1;c_us_2;c_us_3;c_us_4;c_ls_1;c_ls_2;c_ls_3;c_ls_4];
    
    constr_List = [ sos(p_W_low);sos(p_W_up);
        sos(p_CCM);
        sos(su_1); sos(su_2);
        sos(su_3); sos(su_4);
        sos(sl_1); sos(sl_2);
        sos(sl_3); sos(sl_4);        
        W_bounds;
        LQR_constr_W;LQR_constr_rho];
    
    options = sdpsettings('solver','mosek','verbose',1,'debug',1);
    SOS_soln = solvesos(constr_List,w_upper-w_lower,options,coeff_List);
    if (SOS_soln.problem == 0)
        %         bounds(ll) = (1/lambda)*sqrt(value(w_upper/w_lower));
        disp('feasible');
        cond_u = cond;
        cond_best = cond;
    else
        cond_l = cond;
    end
end
pause;
clc;
% % 
% Display

% coeff_sol = clean(double(coeff_List),1e-3);
% W_sol = replace(W,coeff_List,coeff_sol);
% rho_sol = replace(rho,coeff_List,coeff_sol);
% % 
% disp('W')
% sdisplay(W_sol)
% 
% disp('rho')
% sdisplay(rho_sol)

disp('Bound achieved');
disp((1/lambda)*sqrt(value(cond_best)));

%% Figure

% Convert into functions
W_cell = cell(3,3);
for i = 1:3
    for j = 1:3
        W_ij = sdisplay(W_sol(i,j));
        W_ij = strcat('@(x1,x2)',W_ij{1});
        W_cell{i,j} = str2func(W_ij);
    end
end
W_mat = @(x1,x2) [W_cell{1,1}(x1,x2), W_cell{1,2}(x1,x2), W_cell{1,3}(x1,x2);
                  W_cell{2,1}(x1,x2), W_cell{2,2}(x1,x2), W_cell{2,3}(x1,x2);
                  W_cell{3,1}(x1,x2), W_cell{3,2}(x1,x2), W_cell{3,3}(x1,x2);];

x1_range = -4:0.2:4;
x2_range = -4:0.2:4;
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
xlabel('X'); ylabel('Y');






