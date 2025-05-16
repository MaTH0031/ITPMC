close all
clear all
clc

import casadi.*

N =100; % number of control intervals

c=0.02;
p1=50;p2=50;
K1=10;K2=10;

uM1=120;
uM2=0;
uM3=8;

eta1=10;
epsion1=0;
d1=0;


x1_min = -inf;x1_max = 20;
x2_min =  -inf;x2_max = 20;
x3_min =  -inf;x3_max = inf;

v1_min = -0.1;v1_max = 1;
v2_min = -1;v2_max = 1;
v3_min = -1;v3_max = 1;

u1_min = 0;u1_max = 50;
u2_min = 0;u2_max = 50;


x = SX.sym('x',6); % parameters
xdot = SX.sym('xdot',6); % parameters
u = SX.sym('u',2); % parameters
dt = 0.2;

x0 =[0;0;pi;0;0;0];
xf = [-15;-15;-pi/2;0;0;0];

m11=15.8;
m22=23.8;
m33=0.965;
distance=0.26;
M=[m11 0 0;0 m22 0;0 0 m33];
Minv=M^(-1);


R=[cos(x(3)),-sin(x(3)),0;
   sin(x(3)), cos(x(3)),0;
   0,0,1];
% C=[0 0 -b*x(6)-m*x(5)+yvd*x(5)+yrd*x(6);
%    0 0 m*x(4)-xud*x(4);
%    b*x(6)+m*x(5)-yvd*x(5)-yrd*x(6) -m*x(4)+xud*x(4) 0];
% D=[-xu-sign(x(4))*xuu*x(4) 0 0;
%    0 -yv-sign(x(5))*yvv*x(5)-yrv*abs(x(6)) -yr-yvr*abs(x(5))-sign(x(6))*yrr*x(6);
%    0 -nv-nvv*sign(x(5))*x(5)-nrv*abs(x(6)) -nr-nvr*abs(x(5))-nrr*sign(x(6))*x(6)];
C=[0 0 -m22*x(5);
   0 0 m11*x(4);
   m22*x(5) -m11*x(4) 0];
%C=zeros(3,3);
D=[4.5+sign(x(4))*1.32*x(4)+14.4*x(4)*x(4) 0 0;
   0 13.8+32.3*sign(x(5))*x(5) 0;
   0 0 7.5+9.5*x(6)*x(6)];
%D=eye(3,3);
xdot(1:3) = R*x(4:6);
input=[u(1)+u(2),0,(u(2)-u(1))*distance/2]';
xdot(4:6)= Minv*(-C*x(4:6)-D*x(4:6)+input);

f = Function('f', {x, u}, {xdot});

X = SX.sym('X',length(x),N+1);
U = SX.sym('U',length(u),N);
P = SX.sym('P',length(x)*2);

xk = X(:,1); % initial state
g = xk-P(1:length(x)); % initial condition

W1=diag([5+3*sin(2*xf(3)-pi/2) 5+3*sin(-2*xf(3)+pi/2) 0 0 0 0]);
W2=diag([0.001 0.001]);
W3=diag([10+7*sin(2*xf(3)-pi/2) 10+7*sin(-2*xf(3)+pi/2) 0 0 0 0]);

obj=0;
obs=[-5,-1]';
for k = 1:N
    xk = X(:,k); uk = U(:,k);
    xk_next = X(:,k+1);

    err=xk-xf;
    obj = obj+err'*W1*err+uk'*W2*uk;
    obj = obj+5*(atan2(sin(xk(3,1)-xf(3,1)),cos(xk(3,1)-xf(3,1))))^2;


    % 5-th order RK scheme
    k1 = dt*f(xk,uk);
    k2 = dt*f(xk+4*k1/11,uk);
    k3 = dt*f(xk+(9*k1+11*k2)/50,uk);
    k4 = dt*f(xk+(-11*k2+15*k3)/4,uk);
    k5 = dt*f(xk+((81+9*sqrt(6))*k1+(255-55*sqrt(6))*k3+(24-14*sqrt(6))*k4)/600,uk);
    k6 = dt*f(xk+((81-9*sqrt(6))*k1+(255+55*sqrt(6))*k3+(24+14*sqrt(6))*k4)/600,uk);
    xk_next_RK = xk+(4*k1+(16+sqrt(6))*k5+(16-sqrt(6))*k6)/36;
    g = [g;xk_next-xk_next_RK];
end
for k = 1:N
    xk = X(:,k+1); 
    obs_dis=(xk(1,1)-obs(1,1))^2+(xk(2,1)-obs(2,1))^2;

    g = [g;obs_dis];
end

err=X(:,N+1)-xf;
obj = obj+err'*W3*err;
obj = obj+10*(atan2(sin(xk(3,1)-xf(3,1)),cos(xk(3,1)-xf(3,1))))^2;

OPT_variables = [reshape(X,6*(N+1),1);reshape(U,2*N,1)];

nlp_prob = struct('f',obj,'x',OPT_variables,'g',g,'p',P);

% solver options
opts = struct;
%opts.ipopt.linear_solver = 'ma57';
opts.ipopt.max_iter = 500;
opts.ipopt.tol = 1e-5;
opts.ipopt.nlp_scaling_method = 'gradient-based';
opts.ipopt.mehrotra_algorithm = 'no';
opts.ipopt.mu_strategy = 'monotone';
opts.ipopt.hessian_approximation = 'exact'; %'limited-Memory'
opts.ipopt.warm_start_init_point = 'no';
opts.ipopt.print_level =0;
opts.print_time = 0;

solver = nlpsol('solver','ipopt',nlp_prob,opts);
 
% solver input
args = struct;
args.lbx(1:6:6*(N+1),1) = x1_min;
args.ubx(1:6:6*(N+1),1) = x1_max;
args.lbx(2:6:6*(N+1),1) = x2_min;
args.ubx(2:6:6*(N+1),1) = x2_max;
args.lbx(3:6:6*(N+1),1) = x3_min;
args.ubx(3:6:6*(N+1),1) = x3_max;
args.lbx(4:6:6*(N+1),1) = v1_min;
args.ubx(4:6:6*(N+1),1) = v1_max;
args.lbx(5:6:6*(N+1),1) = v2_min;
args.ubx(5:6:6*(N+1),1) = v2_max;
args.lbx(6:6:6*(N+1),1) = v3_min;
args.ubx(6:6:6*(N+1),1) = v3_max;
args.lbx(6*(N+1)+1:2:6*(N+1)+2*N,1) = u1_min;
args.ubx(6*(N+1)+1:2:6*(N+1)+2*N,1) = u1_max;
args.lbx(6*(N+1)+2:2:6*(N+1)+2*N,1) = u2_min;
args.ubx(6*(N+1)+2:2:6*(N+1)+2*N,1) = u2_max;



args.lbg(1:6*(N+1)) = 0;
args.ubg(1:6*(N+1)) = 0;
args.lbg(6*(N+1)+1:6*(N+1)+N) = 4;
args.ubg(6*(N+1)+1:6*(N+1)+N) = inf;

args.p = [x0;xf];

X0 = zeros(length(x),N+1);
U0 = zeros(length(u),N);
args.x0 = [reshape(X0,6*(N+1),1);reshape(U0,2*N,1)];

% solution
tic
sol = solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
toc

x_sol = full(sol.x);
fmin = full(sol.f);

X_sol = reshape(x_sol(1:6*(N+1)),6,(N+1));
U_sol = reshape(x_sol(6*(N+1)+1:6*(N+1)+2*N),2,N);
uf =U_sol(:,N);
u_sol = [U_sol uf];

Ts_inner=1e-1;

x1r_seq_new = X_sol(1,:)';
x2r_seq_new = X_sol(2,:)';

% figure(1)
% subplot(6,1,1)
% plot(t,X_sol(1,:))
% xlabel('t [s]')
% ylabel('X [m]')
% subplot(6,1,2)
% plot(t,X_sol(2,:))
% xlabel('t [s]')
% ylabel('Y [m/s]')
% subplot(6,1,3)
% plot(t,X_sol(3,:))
% xlabel('t [s]')
% ylabel('theta [rad]')
% subplot(6,1,4)
% plot(t,X_sol(4,:))
% xlabel('t [s]')
% ylabel('u [m/s]')
% subplot(6,1,5)
% plot(t,X_sol(5,:))
% xlabel('t [s]')
% ylabel('v [m/s]')
% subplot(6,1,6)
% plot(t,X_sol(6,:))
% xlabel('t [s]')
% ylabel('r [rad/s]')
% 
figure(2)
subplot(3,1,1)
plot(u_sol(1,:))
xlabel('t [s]')
ylabel('U1')
subplot(3,1,2)
plot(u_sol(2,:))
xlabel('t [s]')
ylabel('U2')
subplot(3,1,3)
plot(X_sol(3,:)/3.14*180)
xlabel('t [s]')
ylabel('theta')

figure(3)
%scatter(X_sol(1,:),X_sol(2,:))
scatter(x1r_seq_new,x2r_seq_new)
axis([-5 15 -5 15])
axis square
