function [sys,x0,str,ts,simStateCompliance] = plant(t,x,u,flag,x_initial)
%VSFUNC Variable step S-function example.
%基于无人船动力学的执行器
%预留干扰和参数不确定性通道
%v1 cxl
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(x_initial);

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
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);
  
  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case 1,
    sys = [];

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
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(x_initial)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array
%
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 6;
sizes.NumOutputs     = 6; %先输出全状态
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;     % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [0 0 -pi/4 0 0 0];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = 0.1;      % variable sample time

% speicfy that the simState for this s-function is same as the default
simStateCompliance = 'DefaultSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

global tz
tz = 0.1;
f = RK_fun(u,x,t);
k1 = tz*f;

f = RK_fun(u,x+k1/2,t);
k2 = tz*f;

f = RK_fun(u,x+k2/2,t);
k3 = tz*f;

f = RK_fun(u,x+k3,t);
k4 = tz*f;

xNext = x + (k1+2*k2+2*k3+k4)/6;

sys = xNext;

if (x(1)-15)^2+(x(2)-15)^2<0.009 || x(1)>15.5 || x(2)>15.5 
    sys = x;
end

function sys=mdlOutputs(t,x,u)

sys = x;

function sys=mdlGetTimeOfNextVarHit(t,x,u)

sys = t + 0.1;

function sys=mdlTerminate(t,x,u)

sys = [];

function f=RK_fun(u,x,t)
R=[cos(x(3)) -sin(x(3)) 0;
   sin(x(3)) cos(x(3)) 0;
   0 0 1];
V=[x(4) x(5) x(6)]';
m11=15.8;
m22=23.8;
m33=0.965;
M=[m11 0 0;0 m22 0;0 0 m33];
D=[5.5+1.52*abs(x(4))+17.4*x(4)*x(4) 0 0;
   0 13.8+32.3*abs(x(5)) 0;
   0 0 6.5+7.5*x(6)*x(6)];
% D=[4.5+sign(x(4))*1.32*x(4)+14.4*x(4)*x(4) 0 0;
%    0 13.8+32.3*sign(x(5))*x(5) 0;
%    0 0 7.5+9.5*x(6)*x(6)];
C=[0 0 -m22*x(5);
   0 0 m11*x(4);
   m22*x(5) -m11*x(4) 0];

U=[u(1),u(2),u(3)]';
f=[R*V;
   M^(-1)*(-D*V-C*V+U)];