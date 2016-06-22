clear all; close all; clc;
% clear(yalmip)
%% YALMIP SOS
% x = sdpvar(1,1);y = sdpvar(1,1);
% p = (1+x)^4 + (1-y)^2;
% F = sos(p);
% % sdpsettings('solver','pensdp');
% solvesos(F);

%%

x1 = sdpvar(1,1); x2 = sdpvar(1,1);
y1 = sdpvar(1,1); y2 = sdpvar(1,1);

y = [y1;y2];
f = [-x2-(3/2)*(x1^2)-(1/2)*(x1^3);
       3*x1-x2];
df = [-3*x1 - (3/2)*(x1^2), -1;
       3, -1];

[M_11,c_11,v_11]= polynomial([x1;x2],4);
[M_12,c_12,v_12]= polynomial([x1;x2],4);
[M_22,c_22,v_22]= polynomial([x1;x2],4);

M = [M_11, M_12;
     M_12, M_22];
 
M_dot = [[jacobian(M_11,x1), jacobian(M_11,x2)]*f,...
         [jacobian(M_12,x1), jacobian(M_12,x2)]*f;...
         [jacobian(M_12,x1), jacobian(M_12,x2)]*f,...
         [jacobian(M_22,x1), jacobian(M_22,x2)]*f];

R = -(df.'*M + M*df + M_dot + 0.75*M);

p_M = y'*(M-0.1*eye(2))*y; 
p_R = y'*R*y;
F = [sos(p_M), sos(p_R)];
solvesos(F,norm(c_11,1)+...
           norm(c_12,1)+...
           norm(c_22,1),[],[c_11;c_12;c_22]);
pause;
clc;
%% Display

c_11 = value(c_11);
c_11 = (abs(c_11)>(1e-6)).*c_11;
v_11 = sdisplay(v_11);
v_11 = sym(v_11);
M_11 = sym(c_11','d')*v_11;

c_12 = value(c_12);
c_12 = (abs(c_12)>(1e-6)).*c_12;
v_12 = sdisplay(v_12);
v_12 = sym(v_12);
M_12 = sym(c_12','d')*v_12;

c_22 = value(c_22);
c_22 = (abs(c_22)>(1e-6)).*c_22;
v_22 = sdisplay(v_22);
v_22 = sym(v_22);
M_22 = sym(c_22','d')*v_22;

disp('Ready to verify with SOSTOOLS');
pause;
%% Verify

syms x1 x2
x = [x1;x2];
f = [-x2-(3/2)*(x1^2)-(1/2)*(x1^3);
       3*x1-x2];
df = [-3*x1 - (3/2)*(x1^2), -1;
       3, -1];

prob = sosprogram(x);

M = [M_11, M_12;
     M_12, M_22];
M_dot = diff(M,x1)*f(1) + diff(M,x2)*f(2);
R = -(df.'*M + M*df + M_dot + 0.75*M);

[Q_M,Z_M,H_M] = findsos(M);
[Q_R,Z_R,H_R] = findsos(R);
