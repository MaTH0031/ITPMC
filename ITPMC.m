function [sys,x0,str,ts,simStateCompliance] = ITPMC(t,x,u,flag,K1,K2,p1,p2,ks1,ks2,c,uM1,uM2,uM3,gamma1,gamma2,d1max,d2max,Ts_inner,N)
%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=[];

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=[];

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=[];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 0; 
sizes.NumDiscStates  = 3+8*(N+1);
sizes.NumOutputs     = 14;
sizes.NumInputs      = 6+8*(N+1);
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1; % at least one sample time is needed

sys = simsizes(sizes);


%thetab_hat0 和初始u_controller
x0  = zeros(3+8*(N+1),1);

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';
end
% end mdlInitializeSizes

%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)
m11=15.8;
m22=23.8;
m33=0.965;
M=[m11 0 0;0 m22 0;0 0 m33];
gamma_d=diag([gamma1,gamma2]);
d_max=[d1max d2max]';
d_min=-d_max;

x1=[u(1) u(2) u(3)]';
v=[u(4) u(5) u(6)]';
x1r_seq=u(1+6:N+1+6);
x2r_seq=u(N+2+6:2*(N+1)+6);
x3r_seq=u(2*(N+1)+1+6:3*(N+1)+6);
x4r_seq=u(3*(N+1)+1+6:4*(N+1)+6);
x5r_seq=u(4*(N+1)+1+6:5*(N+1)+6);
x6r_seq=u(5*(N+1)+1+6:6*(N+1)+6);
u1r_seq=u(1+6*(N+1)+6:N+1+6*(N+1)+6);
u2r_seq=u(N+2+6*(N+1)+6:2*(N+1)+6*(N+1)+6);
eventflag=x(3);
if eventflag==0
    x1d=[x1r_seq(1) x2r_seq(1) x3r_seq(1)]';
    vd=[x4r_seq(1) x5r_seq(1) x6r_seq(1)]';
elseif eventflag<=10
    x1d=[x(4+eventflag) x(4+N+1+eventflag) x(4+2*(N+1)+eventflag)]';
    vd=[x(4+3*(N+1)+eventflag) x(4+4*(N+1)+eventflag) x(4+5*(N+1)+eventflag)]';
end

R=[cos(x1(3)),-sin(x1(3)),0;
   sin(x1(3)), cos(x1(3)),0;
   0,0,1];
wx=[0 -v(3) 0;
    v(3) 0 0;
    0 0 0];
Rr=[cos(x1d(3)),-sin(x1d(3)),0;
   sin(x1d(3)), cos(x1d(3)),0;
   0,0,1];
xd_dot=Rr*vd;
z1=R'*(x1-x1d);
z1_dot=-wx*z1+v-R'*xd_dot;
z11=z1(1);
z12=z1(3)+asin(z1(2)/sqrt(1+z1(1)*z1(1)+z1(2)*z1(2)));

kphi1=1.2;
kphi2=0.5;
alpha1=cos(x1(3))*xd_dot(1)+sin(x1(3))*xd_dot(2)-v(3)*z1(2)-K1*z1(1);
temp1=z1_dot(2)-z1(2)*(z1(1)*z1_dot(1)+z1(2)*z1_dot(2))/sqrt(1+z1(1)*z1(1)+z1(2)*z1(2));
alpha2=kphi1*kphi2/sqrt(1+z1(1)*z1(1)+z1(2)*z1(2)-kphi2*kphi2*z1(2)*z1(2))*temp1+vd(3)-K2*z12;

z21=v(1)-alpha1+K1*z11;
z22=v(3)-alpha2+K2*z12;
z2=[z21,z22]';
V=z11*p1*z11+z12*p2*z12+z21*M(1,1)*z21/2+z22*M(3,3)*z22/2;

d=x(1:2);
d_dot=-gamma_d*z2;
d=d+d_dot*Ts_inner;
for i=1:2
   if d(i) > d_max(i)
       d(i) = d_max(i);
     elseif d(i) < d_min(i)
       d(i) = d_min(i);
     else
       d(i) = d(i);
   end
end
sys(1:2)=d;

if V>0.1 || eventflag>=10 || eventflag==0
    eventflag=1;
    sys(4:3+N+1)=x1r_seq;
    sys(4+N+1:3+2*(N+1))=x2r_seq;
    sys(4+2*(N+1):3+3*(N+1))=x3r_seq;
    sys(4+3*(N+1):3+4*(N+1))=x4r_seq;
    sys(4+4*(N+1):3+5*(N+1))=x5r_seq;
    sys(4+5*(N+1):3+6*(N+1))=x6r_seq;
    sys(4+6*(N+1):3+7*(N+1))=u1r_seq;
    sys(4+7*(N+1):3+8*(N+1))=u2r_seq;
else
    eventflag=eventflag+1;
    sys(4:3+8*(N+1))=x(4:3+8*(N+1));
end
sys(3)=eventflag;
end
% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
d=x(1:2);
x1=[u(1) u(2) u(3)]';
v=[u(4) u(5) u(6)]';
x1r_seq=u(1+6:N+1+6);
x2r_seq=u(N+2+6:2*(N+1)+6);
x3r_seq=u(2*(N+1)+1+6:3*(N+1)+6);
x4r_seq=u(3*(N+1)+1+6:4*(N+1)+6);
x5r_seq=u(4*(N+1)+1+6:5*(N+1)+6);
x6r_seq=u(5*(N+1)+1+6:6*(N+1)+6);
u1r_seq=u(1+6*(N+1)+6:N+1+6*(N+1)+6);
u2r_seq=u(N+2+6*(N+1)+6:2*(N+1)+6*(N+1)+6);
eventflag=x(3);
if eventflag==0
    x1d=[x1r_seq(1) x2r_seq(1) x3r_seq(1)]';
    vd=[x4r_seq(1) x5r_seq(1) x6r_seq(1)]';
elseif eventflag<=10
    x1d=[x(4+eventflag) x(4+N+1+eventflag) x(4+2*(N+1)+eventflag)]';
    vd=[x(4+3*(N+1)+eventflag) x(4+4*(N+1)+eventflag) x(4+5*(N+1)+eventflag)]';
end

R=[cos(x1(3)),-sin(x1(3)),0;
   sin(x1(3)), cos(x1(3)),0;
   0,0,1];
wx=[0 -v(3) 0;
    v(3) 0 0;
    0 0 0];
Rr=[cos(x1d(3)),-sin(x1d(3)),0;
   sin(x1d(3)), cos(x1d(3)),0;
   0,0,1];
xd_dot=Rr*vd;
z1=R'*(x1-x1d);
z1_dot=-wx*z1+v-R'*xd_dot;
z11=z1(1);
z12=z1(3)+asin(z1(2)/sqrt(1+z1(1)*z1(1)+z1(2)*z1(2)));

kphi1=1.2;
kphi2=0.5;
alpha1=cos(x1(3))*xd_dot(1)+sin(x1(3))*xd_dot(2)-v(3)*z1(2)-K1*z1(1);
temp1=z1_dot(2)-z1(2)*(z1(1)*z1_dot(1)+z1(2)*z1_dot(2))/sqrt(1+z1(1)*z1(1)+z1(2)*z1(2));
alpha2=kphi1*kphi2/sqrt(1+z1(1)*z1(1)+z1(2)*z1(2)-kphi2*kphi2*z1(2)*z1(2))*temp1+vd(3)-K2*z12;

m11=15.8;
m22=23.8;
m33=0.965;
distance=0.26;
M=[m11 0 0;0 m22 0;0 0 m33];
z21=v(1)-alpha1+K1*z11;
z22=v(3)-alpha2+K2*z12;
V=z11*p1*z11+z12*p2*z12+z21*M(1,1)*z21/2+z22*M(3,3)*z22/2;

us1=-(2*p1-M(1,1)*K1*K1)*z11-(ks1+M(1,1)*K1)*z21;
us2=-(2*p2-M(3,3)*K2*K2)*z12-(ks2+M(3,3)*K2)*z22;

flag_us=1;
flag_d=1;
if eventflag==0 || eventflag==1
    ur_input=[u1r_seq(1)+u2r_seq(1),0,(u2r_seq(1)-u1r_seq(1))*distance/2]';
else
    ur_input=[x(4+6*(N+1)+eventflag)+x(4+7*(N+1)+eventflag),0,(x(4+7*(N+1)+eventflag)-x(4+6*(N+1)+eventflag))*distance/2]';
end
tau_o1=ur_input(1)+flag_us*us1+flag_d*d(1);
tau_o2=ur_input(3)+flag_us*us2+flag_d*d(2);

sys=[tau_o1 0 tau_o2 z11 z12 z21 z22 z1(3) V us1 us2 d(1) d(2) eventflag];

end
% end mdlOutputs
end

