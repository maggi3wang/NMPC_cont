clear all; close all; clc;

%% Constants
n = 2;
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

W = @(x) [0.2543177058, -0.5586044142;
          -0.5586044142,2.492422872-0.006521893711418*x(1)^4];

dW = @(x) {[0, 0;
           0, -0.026087574845672*x(1)^3],...
           zeros(n)};
[geo_Prob,geo_Ke,geo_we,T_e,T_dot_e,geo_Aeq] = ...
        setup_geodesic_calc(n,geodesic_N,W,dW);

%% Setup Auxiliary controller
Q = kron(geo_we',eye(m)); 
Q_bar = Q'*Q;
aux_Prob = setup_opt_aux(geo_Ke,Q_bar,m);

%% Test Auxiliary control computation
start_p = [0;0];
end_p = [-1;0.5];

% Test TOMLAB
tic
[X, X_dot,J_opt,exitflag] = compute_geodesic_tom(geo_Prob,...
    n,geodesic_N,start_p,end_p,T_e,T_dot_e,geo_Aeq);
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
figure(1)
% plot3(start_p(1), start_p(2),start_p(3),'go','markersize',10,'markerfacecolor','g');
plot(start_p(1), start_p(2),'go','markersize',10,'markerfacecolor','g');
hold on
% plot3(end_p(1), end_p(2),end_p(3),'ro','markersize',10,'markerfacecolor','r');
plot(end_p(1), end_p(2),'ro','markersize',10,'markerfacecolor','r');
for k = 1:geo_Ke+1
%     plot3(X(1,k), X(2,k),X(3,k),'ko','markersize',10,'markerfacecolor','k');
    plot(X(1,k), X(2,k),'ko','markersize',10,'markerfacecolor','k');
end
% plot3(X(1,:),X(2,:),X(3,:),'r-','linewidth',2);
plot(X(1,:),X(2,:),'r-','linewidth',2);
hold off
grid on
axis equal
xlabel('x'); ylabel('y');


%Numerics for optimal differential controller
% f = @(x) [-x(1)+x(3); x(1)^2 - x(2) - 2*x(1)*x(3) + x(3);-x(2)];
% df = @(x) [-1, 0, 1;2*x(1)-2*x(3), -1, -2*x(1)+1;0,-1,0];
% B = [0;0;1];

% f = @(x) [-1*x(1) + 2*x(2); 
%           -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];
% df = @(x) [-1, 2;
%            -3, 4-0.75*(x(2)^2)];
% B = [0.5;-2];

f = @(x)[x(2);
        -x(1)+x(2)*sinh(x(1)^2 + x(2)^2)];
df = @(x)[0, 1;
         -1+2*x(1)*x(2)*cosh(x(1)^2 + x(2)^2),sinh(x(1)^2 + x(2)^2)+2*(x(2)^2)*cosh(x(1)^2 + x(2)^2)];
B = [0;
     1];
B_w = [0;0]; w_max = 0;
lambda =  2;


tic
[ctrl_opt,e_rate,exitflag] = compute_opt_aux(aux_Prob,Q,geo_we,geo_Ke,J_opt,X,X_dot,W,dW,...
        f,df,B,zeros(m,1),lambda,n,m);
disp('Got opt control'); disp(ctrl_opt);
toc

%% Simulate non-linear system

ode_options = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);

dt = 0.002;
solve_t = (0:dt:3)';
T_steps = length(solve_t);

t_opt = cell(T_steps-1,1);
% t_s = cell(T_steps-1,1);

state_0_opt = end_p;
state_0_s = end_p;

state_opt = cell(T_steps-1,1);
% state_s = cell(T_steps-1,1);

ctrl_opt = zeros(T_steps,m);
ctrl_s = zeros(T_steps,m);

J_opt = zeros(T_steps,1);
J_s = zeros(T_steps,1);

ctrl_solve_time = zeros(T_steps,2);

solved = ones(T_steps,2);

% rho = @(x) 51.88116196-0.001458799296458*x(1)-...
%        0.001463848144924*x(2)+2.241303600302532*x(2)^2+...
%        7.261640169359435*x(1)^2+1.585568386570665*x(1)*x(2)+...
%        2.655722707839840*x(1)^4+0.611430475833738*x(1)^3*x(2)+...
%        4.105229604485951*x(1)^2*x(2)^2+...
%        2.152384589971836*x(1)*x(2)^3+...
%        6.938732192947612*x(2)^4;

rho =@(x) 5825.652472-1.994747096772980*x(1)-0.169841718823589*x(2)+...
          3.817406236202479e+03*x(1)^2+4.225004246942979e+03*x(2)^2+...
          7.486700415836236e+03*x(1)^4+6.462177245609084e+03*x(1)^2*x(2)^2+...
          9.390587718981402e+03*x(2)^4+18.253220705351527*x(1)^2*x(2)+...
          17.584174648858632*x(2)^3-4.746329713593079e+02*x(1)*x(2)-...
          21.122699043530364*x(1)^3-32.298861002297762*x(1)*x(2)^2-...
          1.930025253840884e+03*x(1)^3*x(2)-2.349756284415636e+03*x(1)*x(2)^3;

cont_sim = 0;
% if (cont_sim)
% [t_s,state_s] = ode113(@(t_s,state_s)ode_sim_c(t_s,state_s,f,B,...
% geo_Prob,W,m,n,geodesic_N,geo_Ke,geo_we,T_e,T_dot_e,geo_Aeq,rho),...
%             [solve_t(1), solve_t(end)],state_0_s,ode_options);
% end
 %% 
for i = 1:T_steps-1
    %Optimal Control
    tic
    [X, X_dot,J_opt(i),solved(i,1)] = compute_geodesic_tom(geo_Prob,...
    n,geodesic_N,[0;0],state_0_opt,T_e,T_dot_e,geo_Aeq);
    ctrl_solve_time(i,1) = toc;
    tic
    [ctrl_opt(i,:),e_rate,solved(i,2)] = compute_opt_aux(aux_Prob,Q,geo_we,geo_Ke,J_opt(i),X,X_dot,W,dW,...
        f,df,B,zeros(m,1),lambda,n,m);
    ctrl_solve_time(i,2) = toc;
    
%     %Standard control
%     [X_s, X_dot_s,J_s(i),exitflag] = compute_geodesic_tom(geo_Prob,...
%     n,geodesic_N,[0;0],state_0_s,T_e,T_dot_e,geo_Aeq);
%     ctrl_s(i,:) = aux_control(rho,X_s,X_dot_s,geo_Ke,geo_we,...
%         B,W,m);
    
    %Simulate Optimal
    [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,ctrl_opt(i,:)',f,B,B_w,w_max),...
        [solve_t(i),solve_t(i+1)],state_0_opt,ode_options);
    t_opt{i} = d_t;
    state_opt{i} = d_state;
    state_0_opt = d_state(end,:)';
    
%     %Simulate Standard
%     [d_t_s,d_state_s] = ode113(@(t,d_state_s)ode_sim(t,d_state_s,ctrl_s(i,:)',f,B,B_w,w_max),...
%         [solve_t(i),solve_t(i+1)],state_0_s,ode_options);
%     t_s{i} = d_t_s;
%     state_s{i} = d_state_s;
%     state_0_s = d_state_s(end,:)';
end

%% If continuous sim - get extra stuff
if (cont_sim)
    ctrl_s = zeros(length(t_s),m);
    J_s = zeros(length(t_s),1);
    for i = 1:length(t_s)
        [X, X_dot,J_s(i),exitflag] = compute_geodesic_tom...
        (geo_Prob,n,geodesic_N,zeros(n,1),state_s(i,:)',...
         T_e,T_dot_e,geo_Aeq);

        ctrl_s(i) = aux_control(rho,X,X_dot,geo_Ke,geo_we,...
            B,W,m);
    end
end

%% plot
close all;
figure()
subplot(2,1,1)
hold on
for i = 1:T_steps-1
    plot(t_opt{i},state_opt{i}(:,1),'r-','linewidth',2);
    plot(t_opt{i},state_opt{i}(:,2),'b-','linewidth',2);
%     plot(t_opt{i},state_opt{i}(:,3),'k-','linewidth',2);
end
grid on; title('Optimized Control');
xlabel('Time [s]'); ylabel('States');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

if (cont_sim)
subplot(2,1,2)
hold on
plot(t_s,state_s,'linewidth',2);
grid on; title('Standard Control');
xlabel('Time [s]'); ylabel('States');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
end

figure()
if (cont_sim)
    stairs(t_s,ctrl_s,'r-','linewidth',2); hold on
else
    stairs(solve_t,ctrl_s,'r-','linewidth',2); hold on
end
stairs(solve_t,ctrl_opt,'b-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); legend('\rho-multiplier','Optimized control');
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
if (cont_sim)
    plot(t_s,J_s,'r-','linewidth',2); hold on
else
    plot(solve_t,J_s,'r-','linewidth',2); hold on
end
plot(solve_t,J_opt(1)*exp(-2*lambda*solve_t),'k--','linewidth',4); 
plot(solve_t,J_opt,'b-','linewidth',2);
grid on
xlabel('Time [s]');
ylabel('Geodesic length');
xlabel('Time [s]'); 
legend('\rho-multiplier','Contraction bound','Optimized control');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
