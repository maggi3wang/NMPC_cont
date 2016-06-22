clear all; close all; clc;

%% Constants
n = 2;
m = 1;

%% Setup NMPC Problem

f = @(x) [-1*x(1) + 2*x(2); 
          -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];
df = @(x) [-1, 2;
           -3, 4-0.75*(x(2)^2)];
B = [0.5;-2];

B_w = [0;1];
w_max = 0.1;

N_mpc = 50;

alpha = 10;
P = [7.9997, -12.2019;
    -12.2019, 27.0777];

Tp = 1.5;
delta = 0.1;
dt = 0.02;

[NMPC_Prob, NMPC_Aeq,L_e,L_e_full] = setup_NMPC(n,m,P,f,...
                    B,df,N_mpc,alpha,Tp,delta,dt);

%% Test MPC Solve

test_state = [3.5;-2.5];
tic
[NMPC_state,NMPC_ctrl,converged] = compute_NMPC(NMPC_Prob,NMPC_Aeq,...
    test_state,n,m,N_mpc,L_e_full);
toc
disp(converged);


% Visualize
close all
figure(); 
plot(NMPC_state(1:(delta/dt),1),...
     NMPC_state(1:(delta/dt),2),'r-','linewidth',2);
hold on
plot(NMPC_state(:,1),NMPC_state(:,2),'b-','linewidth',1);
Ellipse_plot((1/alpha)*P,[0;0],25,'r');
grid on
axis equal
xlabel('X'); ylabel('Y');

pause;

%% Auxiliary controller

K_aux = [-1.3693, 5.1273];

%% Set up non-linear sim


ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

solve_t = (0:dt:10)';
T_steps = length(solve_t);

dt_MPC = delta;
solve_MPC = (0:dt_MPC:10)';
T_steps_MPC = length(solve_MPC);

t = cell(T_steps-1,1);
state_0 = test_state;
state_0_MPC = state_0;

MPC_state = cell(T_steps_MPC-1,1);
MPC_ctrl = cell(T_steps_MPC-1,1);

True_state = cell(T_steps-1,1);
True_ctrl = zeros(T_steps-1,m);

opt_solved = NaN(T_steps-1,1);
ctrl_solve_time = NaN(T_steps-1,1);

i_mpc = 0;

%% Simulate
for i = 1:T_steps-1
    %First Solve MPC
    if (mod(solve_t(i),delta)==0)
        tic
        [NMPC_state,NMPC_ctrl,opt_solved(i,1)] = ...
         compute_NMPC(NMPC_Prob,NMPC_Aeq,state_0_MPC,...
                            n,m,N_mpc,L_e);
        ctrl_solve_time(i,1) = toc;
        
        i_mpc = i_mpc + 1;
        
        MPC_state{i_mpc} = NMPC_state;
        MPC_ctrl{i_mpc} = NMPC_ctrl;
        
        x_nom = MPC_state{i_mpc}(1,:);
        u_nom = MPC_ctrl{i_mpc}(1,:);
        
        %update starting state for next MPC problem
        state_0_MPC = MPC_state{i_mpc}(end,:)';
    else
        i_mpc_use = round((mod(solve_t(i),delta))/dt)+1;
        x_nom = MPC_state{i_mpc}(i_mpc_use,:);
        u_nom = MPC_ctrl{i_mpc}(i_mpc_use,:);
%           x_nom = MPC_state{i_mpc}(1,:);
%           u_nom = MPC_ctrl{i_mpc}(1,:);
    end
    
       
    True_ctrl(i,:) = u_nom+K_aux*(state_0-x_nom');

    %Simulate Optimal
    [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,True_ctrl(i,:)',...
        f,B,B_w,w_max),[solve_t(i),solve_t(i+1)],state_0,ode_options);
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
    plot(t{i},True_state{i}(:,1),'r-','linewidth',2);
    plot(t{i},True_state{i}(:,2),'b-','linewidth',2);
end
grid on
xlabel('Time [s]'); ylabel('States'); legend('x_1','x_2');
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
P_inv = diag([39.0251,486.0402]);
figure()
for i_mpc = 1:T_steps_MPC-1
    plot(MPC_state{i_mpc}(:,1),...
         MPC_state{i_mpc}(:,2),'r-','linewidth',2);
      hold on;
    Ellipse_plot(P_inv,MPC_state{i_mpc}(1,1:2),25,'k');
end
for i = 1:T_steps-1
    plot(True_state{i}(:,1),...
          True_state{i}(:,2),'b-','linewidth',2);    
end
Ellipse_plot((1/alpha)*P,[0;0],25,'r');
xlabel('x_1'); ylabel('x_2');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 


figure()
plot(solve_t(1:end-1),ctrl_solve_time(:,1),'ro','markersize',10,'markerfacecolor','g');
grid on
legend('NMPC');
xlabel('Time [s]');
ylabel('Solve time [s]'); title('Solve time');
xlabel('Time [s]'); 

figure()
plot(solve_t(1:end-1),opt_solved(:,1),'ro','markersize',10,'markerfacecolor','g');
grid on
legend('NMPC');
xlabel('Time [s]');
title('Convergence');

  