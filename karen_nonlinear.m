clear all; close all; clc;

%% Constants
n = 4;
m = 1;
geodesic_N = 4;
%% Setup Geodesic Numerics

%Metric
% W0 = [2.285714286, 0.1428571429, -0.4285714286;
%     0.1428571429, 1.571428571, 0.2857142857;
%     -0.4285714286, 0.2857142857, 1.142857143];
% W1 = [0, -4.571380191053947, 0;
%     -4.571380191053947, -0.571432373322492, 0.857116424029969;
%     0, 0.857116424029969, 0];
% 
% W2 = diag([0;9.142662704421189;0]);
% 
% W = @(x) W0 + W1*x(1) + W2*(x(1)^2);
% dW = @(x) {W1 + 2*W2*x(1),zeros(n,n),zeros(n,n)};

% W = @(x) [2.095084164, -0.4606622251;
%           -0.4606622251, 1.858333021];
% dW = @(x) {zeros(n),zeros(n)};

% W = @(x) [0.2543177058, -0.5586044142;
%           -0.5586044142,2.492422872-0.006521893711418*x(1)^4];
% dW = @(x) {[0, 0;
%            0, -0.026087574845672*x(1)^3],...
%            zeros(n)};
% 
load 'W_nc.mat';
W = @(x) W_0 + W_c*cos(x(1));
dW = @(x) {-W_c*sin(x(1)), zeros(n), zeros(n), zeros(n)};
% dW = @(x) {zeros(n), zeros(n), zeros(n), zeros(n)};


[geo_Prob,geo_Ke,geo_we,T_e,T_dot_e,geo_Aeq] = ...
        setup_geodesic_calc(n,geodesic_N,W,dW);

%% Setup Auxiliary controller
Q = kron(geo_we',eye(m)); 
Q_bar = Q'*Q;
aux_Prob = setup_opt_aux(geo_Ke,Q_bar,m);

%% Test Auxiliary control computation
xstar_init = [0;0;0.15;0];
zstar_init = [0;0;0.15;0];
x_init = [0;0;-0.14;0];
z_init = [0;0;-0.14;0];

% Test TOMLAB
tic
[X, X_dot,J_opt,exitflag] = compute_geodesic_tom(geo_Prob,...
    n,geodesic_N,zstar_init,z_init,T_e,T_dot_e,geo_Aeq);
toc
% Validate solution
Eps =0;
for k = 1:geo_Ke+1
    x_k = X(:,k);
    W_geo = W(x_k);
    M_xdot = W_geo\X_dot(:,k);
    e_k =  (X_dot(:,k)'*M_xdot);
    Eps = Eps + (1/2)*geo_we(k)*(abs(e_k-J_opt)/J_opt);
end

disp('Exit flag'); disp(exitflag);
disp('Eps'); disp(Eps);


%Numerics for optimal differential controller
% f = @(x) [-x(1)+x(3); x(1)^2 - x(2) - 2*x(1)*x(3) + x(3);-x(2)];
% df = @(x) [-1, 0, 1;2*x(1)-2*x(3), -1, -2*x(1)+1;0,-1,0];
% B = [0;0;1];

% f = @(x) [-1*x(1) + 2*x(2); 
%           -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];
% df = @(x) [-1, 2;
%            -3, 4-0.75*(x(2)^2)];
% B = [0.5;-2];

% f = @(x)[x(2);
%         -x(1)+x(2)*sinh(x(1)^2 + x(2)^2)];
% df = @(x)[0, 1;
%          -1+2*x(1)*x(2)*cosh(x(1)^2 + x(2)^2),sinh(x(1)^2 + x(2)^2)+2*(x(2)^2)*cosh(x(1)^2 + x(2)^2)];
% B = [0;
%      1];
% B_w = [0;0]; w_max = 0;
% lambda =  2;

f_aux  = @(z) [z(2)*cos(z(1));
              (z(2)^2)*sin(z(1)) + tan(z(1));
              z(4);
              0];
df_aux = @(z) [-z(2)*sin(z(1)), cos(z(1)), 0, 0;
               -(z(2)^2)*cos(z(1))+(sec(z(1)))^2, 2*z(2)*sin(z(1)), 0, 0;
               0, 0, 0, 1;
               0, 0, 0, 0];
B_aux = [0;-1;0;1];

f= @(x) [x(2);
          sin(x(1));
          x(4);
          0];
B = @(x) [0;
          -cos(x(1));
          0;
          1];
B_w = [0;0;0;1];
      
% lambda  = [0.000236880885503890];
lambda = 0.14;
w_max = 0;       

tic
[ctrl_opt,e_rate,exitflag] = compute_opt_aux(aux_Prob,Q,geo_we,geo_Ke,J_opt,X,X_dot,W,dW,...
        f_aux,df_aux,B_aux,zeros(m,1),lambda,n,m);
disp('Got opt control'); disp(ctrl_opt);
toc

%% Simulate non-linear system

ode_options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

dt = 0.01;
solve_t = (0:dt:20.0)';
T_steps = length(solve_t);

t_opt = cell(T_steps-1,1);

x_0 = x_init;
z_0 = z_init;

state_opt = cell(T_steps-1,1);
ctrl_opt = zeros(T_steps,m);

geod = cell(T_steps-1,2);

J_opt = zeros(T_steps,1);

ctrl_solve_time = zeros(T_steps,2);

solved = ones(T_steps,2);

 %% 
for i = 1:T_steps-1
    %Optimal Control
    disp(solve_t(i));
    tic
    [X, X_dot,J_opt(i),solved(i,1)] = compute_geodesic_tom(geo_Prob,...
    n,geodesic_N,zstar_init,z_0,T_e,T_dot_e,geo_Aeq);
    ctrl_solve_time(i,1) = toc;
    geod{i,1} = X;
    geod{i,2} = X_dot;
    tic
    [ctrl_opt(i,:),e_rate,solved(i,2)] = compute_opt_aux(aux_Prob,Q,geo_we,geo_Ke,J_opt(i),X,X_dot,W,dW,...
        f_aux,df_aux,B_aux,zeros(m,1),lambda,n,m);
    ctrl_solve_time(i,2) = toc;
    
    %Simulate Optimal
    [d_t,d_state] = ode45(@(t,d_state)ode_sim(t,d_state,ctrl_opt(i,:)',f,B,B_w,0),...
        [solve_t(i),solve_t(i+1)],x_0,ode_options);
    t_opt{i} = d_t;
    state_opt{i} = d_state;
    x_0 = d_state(end,:)';
    z_0 = [x_0(1);
           sec(x_0(1))*x_0(2);
           x_0(3);
           x_0(4)];
    
end

%% plot
close all;
figure()
hold on
for i = 1:50:(T_steps-1)
    plot(t_opt{i},state_opt{i}(end,1),'rx','linewidth',2);
    plot(t_opt{i},state_opt{i}(end,2),'bx','linewidth',2);
    plot(t_opt{i},state_opt{i}(end,3),'kx','linewidth',2);
    plot(t_opt{i},state_opt{i}(end,4),'gx','linewidth',2);
    
end
hold off
grid on; title('Optimized Control');
xlabel('Time [s]'); ylabel('States');
legend('\theta','\omega','\xi','v');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)


figure()
stairs(solve_t,ctrl_opt,'b-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

figure()
plot(solve_t,ctrl_solve_time(:,1),'go','markersize',10,'markerfacecolor','g');
hold on
plot(solve_t,ctrl_solve_time(:,2),'ro','markersize',10,'markerfacecolor','r');
grid on
legend('Geodesic computation time','Aux computation');
xlabel('Time [s]');
ylabel('Solve time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

figure()
plot(solve_t,solved(:,1),'go','markersize',10,'markerfacecolor','g');
hold on
plot(solve_t,solved(:,2),'ro','markersize',10,'markerfacecolor','r');
grid on
legend('Geodesic Solved','Aux Solved');
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)


figure()
hold on
plot(solve_t,J_opt(1)*exp(-2*lambda*solve_t),'k--','linewidth',4); 
plot(solve_t,J_opt,'b-','linewidth',2);
hold off
grid on
xlabel('Time [s]');
ylabel('Geodesic energy');
xlabel('Time [s]'); 
legend('Contraction bound','Optimized control');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)