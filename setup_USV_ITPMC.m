clc;
% clear
import casadi.*

%% Hyperparameter definition
N = 20;
Ts_inner=0.1;
p1=10;p2=5;
K1=2;K2=0.7;
ks1=20;ks2=9.5;
gamma1=10;gamma2=10;
d1max=10;d2max=0.5;
uM1=100;uM2=0;uM3=6.6;
epsion11=3;epsion12=0.3;
epsion21=5;epsion22=0.5;
epsion31=0.5;epsion32=0.3;
c=0.1;

W1=diag([8 2 0 0 0 0]);
W2=0.1;
W3=0.001;
W4=diag([0.001 0 0.001]);
W5=diag([17 3 0 0 0 0]);


%% Constraints definition
x1_min = -inf;x1_max = 15.1;
x2_min =  -inf;x2_max = 15.1;
x3_min =  -inf;x3_max = inf;

v1_min = -0.1;v1_max = 1.5-epsion31;
v2_min = -1;v2_max = 1;
v3_min = -1+epsion32;v3_max = 1-epsion32;

u1_min = 0.1;u1_max = (uM1-epsion11-epsion21-3*d1max)/2;
u2_min = u1_min;u2_max = u1_max;

%% Casadi variables
x = SX.sym('x',6);
xdot = SX.sym('xdot',6);
u = SX.sym('u',2);
dt=0.1;

x0 = MX.sym('x0',6);
xf = MX.sym('xf',6);

m11=15.8;
m22=23.8;
m33=0.965;
distance=0.26;
M=[m11 0 0;0 m22 0;0 0 m33];
Minv=M^(-1);

R=[cos(x(3)),-sin(x(3)),0;
   sin(x(3)), cos(x(3)),0;
   0,0,1];
C=[0 0 -m22*x(5);
   0 0 m11*x(4);
   m22*x(5) -m11*x(4) 0];
D=[4.5+sign(x(4))*1.32*x(4)+14.4*x(4)*x(4) 0 0;
   0 13.8+32.3*sign(x(5))*x(5) 0;
   0 0 7.5+9.5*x(6)*x(6)];
xdot(1:3) = R*x(4:6);
input=[u(1)+u(2),0,(u(2)-u(1))*distance/2]';
xdot(4:6)= Minv*(-C*x(4:6)-D*x(4:6)+input);

f = Function('f', {x, u}, {xdot});
X = SX.sym('X',length(x),N+1);
U = SX.sym('U',length(u),N);
P = SX.sym('P',length(x)*2);

%% Construct optimization problem 
obj=0;

xk = X(:,1);
g = xk-P(1:length(x));

for k = 1:N
    xk = X(:,k); uk = U(:,k);
    xk_next = X(:,k+1);
    
    err=xk-P(length(x)+1:2*length(x));
    obj = obj+err'*W1*err;
    obj = obj+5*(atan2(sin(xk(3,1)-P(length(x)+3)),cos(xk(3,1)-P(length(x)+3))))^2;
    obj = obj+W2*(v1_max-xk(4,1));
    obj = obj+W3*(uk(1,1)^1.5+uk(2,1)^1.5);
    inputhere=[uk(1,1)+uk(2,1),0,(uk(2,1)-uk(1,1))*distance/2]';
    obj = obj+inputhere'*W4*inputhere;
    
    k1 = dt*f(xk,uk);
    k2 = dt*f(xk+k1/2,uk);
    k3 = dt*f(xk+k2/2,uk);
    k4 = dt*f(xk+k3,uk);
    xk_next_RK = xk+(k1+2*k2+2*k3+k4)/6;
    
    g = [g;xk_next-xk_next_RK];
end
err=xk_next-P(length(x)+1:2*length(x));
obj = obj+err'*W5*err;
obj = obj+10*(atan2(sin(xk(3,1)-P(length(x)+3)),cos(xk(3,1)-P(length(x)+3))))^2;

%% solver options
OPT_variables = [reshape(X,6*(N+1),1);reshape(U,2*N,1)];
nlp_prob = struct('f',obj,'x',OPT_variables,'g',g,'p',P);
opts = struct;
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
 
%% solver input
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
args.p = [x0;xf];
X0 = zeros(length(x),N+1);
U0 = zeros(length(u),N);
args.x0 = [reshape(X0,6*(N+1),1);reshape(U0,2*N,1)];

%% obtain solution
sol = solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
x_sol = full(sol.x);
fmin = full(sol.f);
X_sol = reshape(x_sol(1:6*(N+1)),6,(N+1));
U_sol = reshape(x_sol(6*(N+1)+1:6*(N+1)+2*N),2,N);
x1_sol = X_sol(1,:)';x2_sol = X_sol(2,:)';x3_sol = X_sol(3,:)';x4_sol = X_sol(4,:)';x5_sol = X_sol(5,:)';x6_sol = X_sol(6,:)';
u1_sol = U_sol(1,:)';u2_sol = U_sol(2,:)';
u_sol = [u1_sol;
         u1_sol(end);
         u2_sol;
         u2_sol(end)];

%% establish dynamic links
Planning = Function('Planning',{x0,xf},{x1_sol,x2_sol,x3_sol,x4_sol,x5_sol,x6_sol,u_sol});
Planning.save('Planning.casadi');
lib_path = GlobalOptions.getCasadiPath();
inc_path = GlobalOptions.getCasadiIncludePath();
mex('-v',['-I' inc_path],['-L' lib_path],'-lcasadi', 'Planning.c')
