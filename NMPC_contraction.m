clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants
n = 4;
m = 1;

%% Setup Geodesic Numerics

lambda =  0.25;
load 'W_nc.mat';
W = @(x) W_0 + W_c*cos(x(1));
dW = @(x) {-W_c*sin(x(1)), zeros(n), zeros(n), zeros(n)};

geodesic_N = 4;
[geo_Prob,geo_Ke,geo_we,T_e,T_dot_e,geo_Aeq] = ...
        setup_geodesic_calc(n,geodesic_N,W,dW);

%Assemble geodesic struct for MPC
geodesic_MPC = struct('geo_Prob',geo_Prob,'geodesic_N',geodesic_N,'T_e',T_e,'T_dot_e',T_dot_e,...
                      'geo_Aeq',geo_Aeq);

%% Setup NMPC Problem

f_aux  = @(z) [z(2)*cos(z(1));
              (z(2)^2)*sin(z(1)) + tan(z(1));
              z(4);
              0];
df_aux = @(z) [-z(2)*sin(z(1)), cos(z(1)), 0, 0;
               -(z(2)^2)*cos(z(1))+(sec(z(1)))^2, 2*z(2)*sin(z(1)), 0, 0;
               0, 0, 0, 1;
               0, 0, 0, 0];
B_aux = [0;-1;0;1];

B_w = 0.5*B_aux;

w_max = 0.002;
RPI_bound = (w_max*40)^2;

N_mpc = 180;

A_LQR = [0, 1, 0, 0;
         1, 0, 0, 0;
         0, 0, 0, 1;
         0, 0, 0, 0];
B_LQR = B_aux;
Q_LQR = diag([1,1,15,1]);
R_LQR = 10;

[K_LQR, P, E_LQR] = lqr(A_LQR,B_LQR,Q_LQR,R_LQR);
alpha = min(eig(P))*(0.8^2)/(max(eig(K_LQR'*K_LQR)));


Tp = 18;
delta = 2;
dt = 0.01;

state_constr = [-0.8*pi/4;-0.4;-1;-0.5]; %theta,eta,cart_p,cart_v
ctrl_constr = -1*ones(m,1);

zstar = [0;0;0.2;0];
[NMPC_Prob, NMPC_Aeq,L_e,L_e_full] = setup_NMPC(n,m,f_aux,...
           B_aux,df_aux,state_constr,ctrl_constr,N_mpc,Tp,delta,dt,P,alpha,W,geodesic_MPC,RPI_bound,zstar);

%% Test MPC Solve

test_state = [0;-0.05;-0.25;-0.1];
tic
[NMPC_state,NMPC_ctrl,converged] = compute_NMPC(NMPC_Prob,NMPC_Aeq,...
    test_state,test_state,state_constr,zstar,n,m,N_mpc,L_e_full);
toc
disp(converged);

% pause;
% Visualize
close all
figure(); 
% 
hold on
% ellipsoid_slice(M0,RPI_bound,[3;1],[2;4],test_state([2;4]),0,NMPC_state(1,:)',100,'g');
% ellipsoid_slice(M0,RPI_bound,[3;1],[2;4],test_state([2;4]),1,NMPC_state(1,:)',100,'g');
% 
% ellipsoid_slice(P,alpha,[3;1],[2;4],NMPC_state(end,[2,4])',0,zstar,25,'r');
% ellipsoid_slice(P,alpha,[3;1],[2;4],NMPC_state(end,[2,4])',1,zstar,25,'r');
plot(NMPC_state(1:(delta/dt),3),...
     NMPC_state(1:(delta/dt),1),'r-','linewidth',2);
plot(NMPC_state(:,3),NMPC_state(:,1),'b-','linewidth',1);
plot(test_state(3),test_state(1),'go','markersize',15,'markerfacecolor','g');
plot(NMPC_state(1,3),NMPC_state(1,1),'ro','markersize',15,'markerfacecolor','r');
grid on
axis equal
xlabel('\xi'); ylabel('\theta');

pause;

%% Test Geodesic Numerics

tic
[X, X_dot,J_opt,exitflag] = compute_geodesic_tom(geo_Prob,...
    n,geodesic_N,NMPC_state(1,:)',test_state,T_e,T_dot_e,geo_Aeq);
toc;
disp(exitflag);

                
%% Setup Auxiliary controller
Q = kron(geo_we',eye(m)); 
Q_bar = Q'*Q;
aux_Prob = setup_opt_aux(geo_Ke,Q_bar,m);

tic
[ctrl_opt,e_rate,exitflag] = compute_opt_aux(aux_Prob,Q,geo_we,geo_Ke,J_opt,X,X_dot,W,dW,...
        f_aux,df_aux,B_aux,zeros(m,1),lambda,n,m);
toc;
disp(exitflag);
%% Set up non-linear sim

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

t_end = 25;
solve_t = (0:dt:t_end)';
T_steps = length(solve_t);

dt_MPC = delta;
solve_MPC = (0:dt_MPC:t_end)';
T_steps_MPC = length(solve_MPC);

t = cell(T_steps-1,1);
state_0 = test_state;
state_0_MPC = state_0;

MPC_state = cell(T_steps_MPC-1,1);
MPC_ctrl = cell(T_steps_MPC-1,1);

Geod = cell(T_steps-1,1);
Aux_ctrl = zeros(T_steps-1,m);

True_state = cell(T_steps-1,1);
True_ctrl = zeros(T_steps-1,m);

ctrl_solve_time = zeros(T_steps-1,3);
ctrl_solve_time(:,1) = NaN;

opt_solved = zeros(T_steps-1,3);
opt_solved(:,1) = NaN;

sep_dist = zeros(T_steps-1,2);
sep_dist(:,2) = NaN;

i_mpc = 0;

      
%% Simulate

for i = 1:T_steps-1
    fprintf('%d/%d \n',i,T_steps-1);
    %First Solve MPC
    if (mod(solve_t(i),delta)==0)
%         sep_dist(i,1) = (state_0-state_0_MPC)'*M0*(state_0-state_0_MPC);
        [~, ~,J_opt,~] = compute_geodesic_tom(geo_Prob,...
            n,geodesic_N,state_0_MPC,state_0,T_e,T_dot_e,geo_Aeq);
        sep_dist(i,1) = J_opt;
        
        tic
        [NMPC_state,NMPC_ctrl,opt_solved(i,1)] = ...
         compute_NMPC(NMPC_Prob,NMPC_Aeq,state_0,...
                            state_0_MPC,state_constr,zstar,n,m,N_mpc,L_e);
        ctrl_solve_time(i,1) = toc;
        
        i_mpc = i_mpc + 1;
        
        MPC_state{i_mpc} = NMPC_state;
        MPC_ctrl{i_mpc} = NMPC_ctrl;
        
        x_nom = MPC_state{i_mpc}(1,:);
        u_nom = MPC_ctrl{i_mpc}(1,:);
        
        [~, ~,J_opt,~] = compute_geodesic_tom(geo_Prob,...
            n,geodesic_N,x_nom',state_0,T_e,T_dot_e,geo_Aeq);
        
%         sep_dist(i,2) = (state_0-x_nom')'*M0*(state_0-x_nom');
        sep_dist(i,2) = J_opt;
        
        %update starting state for next MPC problem
        state_0_MPC = MPC_state{i_mpc}(end,:)';
    else
        i_mpc_use = round((mod(solve_t(i),delta))/dt)+1;
        x_nom = MPC_state{i_mpc}(i_mpc_use,:);
        u_nom = MPC_ctrl{i_mpc}(i_mpc_use,:);
    end
    
    %Optimal Control
    tic
    [X, X_dot,J_opt,converged] = compute_geodesic_tom(geo_Prob,...
    n,geodesic_N,x_nom',state_0,T_e,T_dot_e,geo_Aeq);
    ctrl_solve_time(i,2) = toc;
    opt_solved(i,2) = converged;
    Geod{i} = X';
    tic
    [Aux_ctrl(i,:),opt_solved(i,3)] = compute_opt_aux(aux_Prob,Q,geo_we,geo_Ke,J_opt,X,X_dot,W,dW,...
        f_aux,df_aux,B_aux,u_nom,lambda,n,m);
    ctrl_solve_time(i,3) = toc;   

    
    True_ctrl(i,:) = u_nom+Aux_ctrl(i,:);
%     sep_dist(i,1) = (state_0-x_nom')'*M0*(state_0-x_nom');
    sep_dist(i,1) = J_opt;
    

    %Simulate Optimal
%     w_dist = 0.0*w_max + (w_max)*sin(2*pi*0.5*solve_t(i));
    w_dist = w_max;
    [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,True_ctrl(i,:)',...
        f_aux,B_aux,B_w,w_dist),[solve_t(i),solve_t(i+1)],state_0,ode_options);
    t{i} = d_t;
    True_state{i} = d_state;
    state_0 = d_state(end,:)';
    
end

%% Plot
close all

% State Trajectory
figure()
hold on
for i = 1:T_steps-1
    plot(t{i},True_state{i}(:,3),'r-','linewidth',2);
    plot(t{i},True_state{i}(:,1),'b-','linewidth',2);
end
grid on
xlabel('Time [s]'); ylabel('States'); legend('\xi','\theta');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Control Trajectory
figure()
for i_mpc = 1:T_steps_MPC-1
    stairs([(i_mpc-1)*delta:dt:i_mpc*delta]',...
             MPC_ctrl{i_mpc},'b-','linewidth',2); hold on
end
stairs(solve_t(1:end-1),True_ctrl,'r-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

%2D State Plot
figure()
for i_mpc = 1:T_steps_MPC-1
    plot(MPC_state{i_mpc}(:,3),...
          MPC_state{i_mpc}(:,1),'r-','linewidth',2);
    hold on;   
end
for i = 1:T_steps-1
    plot(True_state{i}(:,3),...
          True_state{i}(:,1),'b-','linewidth',2);
%     plot(Geod{i}(:,1),Geod{i}(:,2),'k-','linewidth',2);
end
xlabel('x_1'); ylabel('x_2');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 

figure()
% plot(solve_t(1:end-1),ctrl_solve_time(:,1),'ro','markersize',10,'markerfacecolor','g');
hold on
% plot(solve_t(1:end-1),ctrl_solve_time(:,2),'bs','markersize',10,'markerfacecolor','m');
plot(solve_t(1:end-1),ctrl_solve_time(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
% legend('Aux');
xlabel('Time [s]');
ylabel('Solve time [s]'); %title('Solve time');
xlabel('Time [s]'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

figure()
plot(solve_t(1:end-1),opt_solved(:,1),'ro','markersize',20,'markerfacecolor','g');
hold on
plot(solve_t(1:end-1),opt_solved(:,2),'bs','markersize',15,'markerfacecolor','m');
plot(solve_t(1:end-1),opt_solved(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('NMPC','Geodesic','Aux');
xlabel('Time [s]');
title('Convergence');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
 
figure()
plot(solve_t(1:end-1),sep_dist(:,1),'b-','linewidth',2);
hold on
plot(solve_t(1:end-1),RPI_bound*ones(T_steps-1,1),'r-','linewidth',2);
grid on
legend('Distance between trajectories','RCI bound');
xlabel('Time [s]');
plot(solve_t(1:end-1),sep_dist(1:end,2),'bo','markersize',10,'markerfacecolor','b','linewidth',2);


