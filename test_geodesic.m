clear all; close all; clc;

%% NUMERICS SETUP

%Method
method = 2; %Lagrange method = 1

%Chebyshev polynomial order
N = 6;

%Collocation points
if (method==1)
    K = N;
else
    K = N+2;
end

%state dimension
n = 3;
%control dimension
m = 3;

solver_options = optimset(...
    'Algorithm','sqp',...
    'GradObj','on',...
    'Display', 'off',...
    'TolFun', 1e-6,...
    'TolCon', 1e-6,...
    'TolX', 1e-6);
offlineOptions.solver.MaxIter = 20;


%Obtain Checyschev Pseudospectral Numerics
cheby_numerics = setup_ChebyschevPseudo(N,K,n,m,method);
Phi = zeros(n,n*(N+1),K+1);
Phi_dot = zeros(n,n*(N+1),K+1);
for k = 1:K+1
    Phi(:,:,k) = kron(eye(n),cheby_numerics.T(:,k)');
    Phi_dot(:,:,k) = 2*kron(eye(n),cheby_numerics.T_dot(:,k)');
end


%% Metric
W0 = [2.286 0.143 -0.429;
    0.143 1.571 0.286;
    -0.429 0.286 1.143];
W1 = [0, -4.571, 0;
    -4.571, -0.572, 0.857;
    0, 0.857, 0];

W2 = diag([0;9.143;0]);

%% Cost function

if (method==1)
    geodesic_cost_fnc =  @(vars) Geodesic_cost2(vars,cheby_numerics.w,...
        cheby_numerics.D,N,n,W0,W1,W2);
else
    geodesic_cost_fnc = @(vars) Geodesic_cost(vars,cheby_numerics.w,K,N,n,...
        cheby_numerics.T,...
        W0,W1,W2,Phi,Phi_dot);
end

%% Problem def

% Linear equality constraints
Aeq = cheby_numerics.A;

start_p = [0;0;0];
each_end = ones(1,10);
End_p = [];
for j = 1:10
    End_p = [End_p, kron(each_end,[j;j;j])];
end

run_time = zeros(100,1); Eps = zeros(100,1);
for iter = 1:100
    end_p = End_p(:,iter);
    beq = [start_p; end_p];
    
    % path_guess = pinv(Aeq)*beq;
    vars_0 = zeros(n*(N+1),1);
    for i = 1:n
        vars_0((i-1)*(N+1)+1:(i-1)*(N+1)+2) = [(start_p(i)+end_p(i))/2;
            -(start_p(i)-end_p(i))/2];
    end
    
    geodesic_problem = createOptimProblem('fmincon','x0',vars_0,...
        'objective',geodesic_cost_fnc,'Aeq',Aeq,'beq',beq,'options',solver_options);
    
    %%
    tic
    [vars_opt, J_opt, exitflag, solveroutput] = fmincon(geodesic_problem);
    run_time(iter) = toc;
    
    disp('Exit flag'); disp(exitflag);
    
    if (method==1)
        X = vars_opt;
        X_dot = 2*cheby_numerics.D*X;
        X = reshape(X,n,N+1);
        X_dot = reshape(X_dot,n,N+1);
    else
        C_opt = (reshape(vars_opt,N+1,n))';
        X = C_opt*cheby_numerics.T;
        X_dot = 2*C_opt*cheby_numerics.T_dot;
    end
    
    
    % Validate solution
    for k = 1:K+1
        x_k = X(:,k);
        W = W0 + W1*x_k(1) + W2*(x_k(1)^2);
        %     W = eye(3);
        M_xdot = W\X_dot(:,k);
        e_k =  (X_dot(:,k)'*M_xdot);
        Eps(iter) = Eps(iter) + (1/2)*cheby_numerics.w(k)*(abs(e_k-J_opt)/J_opt);
    end
    
end

%% PLOT PERFORMANCE
figure()
subplot(2,1,1)
plot(1:100,run_time,'b-','linewidth',2);
grid on

subplot(2,1,2)
plot(1:100,Eps,'b-','linewidth',2);
grid on

%% PLOT SOLUTIOn

plot_curve = 0;

if (plot_curve)
    
    disp('Energy'); disp(J_opt);
    disp('Validation Eps'); disp(Eps);
    
    % Plot
    
    figure()
    plot3(start_p(1), start_p(2),start_p(3),'go','markersize',10,'markerfacecolor','g');
    hold on
    plot3(end_p(1), end_p(2),end_p(3),'ro','markersize',10,'markerfacecolor','r');
    for k = 1:K+1
        plot3(X(1,k), X(2,k),X(3,k),'ko','markersize',10,'markerfacecolor','k');
    end
    plot3(X(1,:),X(2,:),X(3,:),'r-','linewidth',2);
    hold off
    grid on
    axis equal
    
end

