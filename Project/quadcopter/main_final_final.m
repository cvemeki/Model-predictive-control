% loads:
%    hovering equilibrium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outer controller optimizer instance
clc
clear all
load('quadData.mat')
outerController = getOuterController(Ac, 'cplex');
disp('Data successfully loaded')
%General configuration
%initialization
d2r = 2*pi/360; %convert deg to rad
T = 10; %simulation time
N=3/sys.Ts; %horizon length
% Define constraints limits
zDotMax = 1;
alphaMax = 10*d2r;
betaMax = 10*d2r;
alphaDotMax = 15*d2r;
betaDotMax = 15*d2r;
gammaDotMax = 60*d2r;
uMax = 1-us(1);
uMin = -us(1);
% H*x<=h
Hx = [-1 0 0 0 0 0 0;1 0 0 0 0 0 0;
      0 -1 0 0 0 0 0;0 1 0 0 0 0 0;
      0 0 -1 0 0 0 0;0 0 1 0 0 0 0;
      0 0 0 0 0 0 0;0 0 0 0 0 0 0;%gamma no constraint
      0 0 0 0 -1 0 0;0 0 0 0 1 0 0;
      0 0 0 0 0 -1 0;0 0 0 0 0 1 0;
      0 0 0 0 0 0 -1;0 0 0 0 0 0 1;];
hx = [zDotMax;zDotMax;alphaMax;alphaMax;betaMax;betaMax;0;0;
      alphaDotMax;alphaDotMax;betaDotMax;betaDotMax;gammaDotMax;gammaDotMax;];
% Hu*u<=hu
Hu = [-1 0 0 0;
      1 0 0 0;
      0 -1 0 0;
      0 1 0 0;
      0 0 -1 0;
      0 0 1 0;
      0 0 0 -1;
      0 0 0 1;];
hu = [-uMin;uMax;-uMin;uMax;-uMin;uMax;-uMin;uMax;];
%definition of the initial state Z_rate alpha beta gamma alpha_rate beta_rate gamma_rate 
x0 = [-1 10*d2r -10*d2r 120*d2r 0 0 0]';
%% %%%%%%%%%%%%%% First MPC controller %%%%%%%%%%%%%%%%%%%
yalmip('clear');
close all;
Q = diag([50, 10000, 10000, 2, 0.1, 0.1, 0.1]);
R = eye(4);
[K,P] = dlqr(sys.A,sys.B,Q,R); %K=-Klqr
%compute the maximum invariant set
X = Polyhedron([Hx;Hu*-K],[hx;hu]);
while(1)
    preX = Polyhedron(X.A*(sys.A-sys.B*K),X.b);
    X_inter = Polyhedron([preX.A;X.A], [preX.b;X.b]);
    if X_inter == X
        X_inf = X_inter;
        break;
    end
    X = X_inter;
end
Xf = X_inf;
%Define optimization variables  
x = sdpvar(7,N,'full');
u = sdpvar(4,N,'full');
%Define constraints and objective 
con = ( x(:,2) == sys.A*x(:,1) + sys.B*u(:,1)) + (Hu*u(:,1) <= hu);
obj = u(:,1)'*R*u(:,1);
for i = 2:N-1
    con = con + (x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i));   %System dynamics
    con = con + (Hx*x(:,i) <= hx);                   %State constraints
    con = con + (Hu*u(:,i) <= hu);                   %Input constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); %Cost function
end
con = con + (Xf.A*x(:,N) <= Xf.b);   %Terminal constraint
obj = obj + x(:,N)'*P*x(:,N);      %Terminal weight 
%Compile the matrices
innerController = optimizer(con,obj,[],x(:,1),u(:,1));
%Can now compute the optimal control input using
simQuad( sys, innerController, x0, T);
%% Reference tracking - no disturbance, no invariant sets
fprintf('PART II - reference tracking...\n')
yalmip('clear');
close all;
Q = diag([50, 10000, 10000, 2, 0.1, 0.1, 0.1]);
R = eye(4);
[K,P] = dlqr(sys.A,sys.B,Q,R); %K=-Klqr
% Compute matrix for steady state
C = [eye(4) zeros(4,3)];
steady_mat = [eye(7)-sys.A -sys.B; C zeros(4,4)];

% Define optimization variables
x = sdpvar(7,N,'full');
u = sdpvar(4,N,'full');
ref = sdpvar(4,1,'full');
steady_x = sdpvar(7,1,'full');
steady_u = sdpvar(4,1,'full');

%Define constraints and objective
con = [];
con = con + (steady_mat*[steady_x;steady_u;] ==  [zeros(7,1); ref]);
con = con + (x(:,2)-steady_x == sys.A*(x(:,1)-steady_x) + sys.B*(u(:,1)-steady_u)); % System dynamics
con = con + (Hu*u(:,1) <= hu); % Input constraints
con = con + (Hx*x(:,1) <= hx); % State constraints
obj=(x(:,1)-steady_x)'*Q*(x(:,1)-steady_x) + (u(:,1)-steady_u)'*R*(u(:,1)-steady_u);
for i = 2:(N-1)
    con = con + (x(:,i+1)-steady_x == sys.A*(x(:,i)-steady_x) + sys.B*(u(:,i)-steady_u)); % System dynamics
    con = con + (Hu*u(:,i) <= hu);% Input constraints
    con = con + (Hx*x(:,i) <= hx);% State constraints
    obj = obj + (x(:,i)-steady_x)'*Q*(x(:,i)-steady_x) + (u(:,i)-steady_u)'*R*(u(:,i)-steady_u); % Cost function
end
obj = obj + (x(:,N)-steady_x)'*P*(x(:,N)-steady_x); % Terminal weight
% Compute optimizer
innerController = optimizer(con, obj, [], [x(:,1);ref], u(:,1));

%Constant reference
reference = [0.4; 2*d2r; 2*d2r; 40*d2r];
% simQuad(sys, innerController, x0, T, reference);
%Varying reference
m=T/sys.Ts;
ref1 = [0.5; 3*d2r; 3*d2r; 30*d2r];
varyRef = 5*[0.05*sin((1:m)/m*2*pi); 0.01*cos((1:m)/m*2*pi); 0.01*sin((1:m)/m*2*pi); 0.1*sin((1:m)/m*2*pi)];
reference_vary = zeros(4, m);
for i=1:m
    reference_vary(:,i) = ref1 + varyRef(:,i);
end
simQuad( sys, innerController, x0, T, reference_vary);
%% Nonlinear model simulation - no disturbance
fprintf('Running the FIRST NL model simulation...\n')
sim('simulation1.mdl') 

% pause
%% Disturbance estimation
%estimator
Bd = eye(7);
A_est = [sys.A Bd; zeros(7) eye(7)];
B_est = [sys.B; zeros(7,4)];
Cd = zeros(4,7);
C_est = [eye(7) zeros(7)];
% Compute matrix for steady state
C = [eye(4) zeros(4,3)];
steady_mat = [eye(7)-sys.A -sys.B; C zeros(4,4)];

pole = linspace(0.88,0.91,14);
L = place(A_est',C_est', pole)'; % Observer
filter.Af = A_est-L*C_est;
filter.Bf = [B_est L];
% Offset free MPC
fprintf('PART III - OFFSET FREE / Disturbance rejection...\n')
yalmip('clear');
close all;
T = 10; %simulation time
N=20;
% Q = diag([250, 1000, 1000, 5, 0.1, 0.1, 1]);% corresponding to first poles
Q = diag([250, 10000, 10000, 5, 0.1, 0.1, 5]);%For big pole values
% Q = diag([5, 200, 200, 10, 10, 10, 1]);
R = eye(4);
[K,P] = dlqr(sys.A,sys.B,Q,R); %K=-Klqr
% Declare variables
x = sdpvar(7,N,'full');
u = sdpvar(4,N,'full');
ref = sdpvar(4,1,'full');
steady_x = sdpvar(7,1,'full');
steady_u = sdpvar(4,1,'full');
dest = sdpvar(7,1,'full');
% Define constraints and objective
con = [];
con = con + (steady_mat*[steady_x;steady_u] ==  [Bd*dest; ref-Cd*dest]);
con = con + (x(:,2) == sys.A*x(:,1) + sys.B*u(:,1) + dest); % System dynamics

con = con + (-zDotMax<=x(1,:)<=zDotMax);
con = con + (-alphaMax<=x(2,:)<=alphaMax); 
con = con + (-betaMax<=x(3,:)<=betaMax);
con = con + (-alphaDotMax<=x(5,:)<=alphaDotMax);
con = con + (-betaDotMax<=x(6,:)<=betaDotMax);
con = con + (-gammaDotMax<=x(7,:)<=gammaDotMax);
con = con + (uMin<=u<=uMax);
obj=(x(:,1)-steady_x)'*Q*(x(:,1)-steady_x) + (u(:,1)-steady_u)'*R*(u(:,1)-steady_u);
for i = 2:(N-1)
    con = con + (x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i) + dest); % System dynamics
    obj = obj + (x(:,i)-steady_x)'*Q*(x(:,i)-steady_x) + (u(:,i)-steady_u)'*R*(u(:,i)-steady_u); % Cost function
end
obj = obj + (x(:,N)-steady_x)'*P*(x(:,N)-steady_x); % Terminal weight

% Compute optimizer
option = sdpsettings('solver','cplex');
innerController = optimizer(con, obj, option, [x(:,1);ref;dest], u(:,1));

%Constant reference
reference = [0.4; 2*d2r; -2*d2r; 60*d2r];

% simQuad(sys, innerController, x0, T, reference,filter);

%Varying reference
m=T/sys.Ts;
varyRef = [0.05*sin((1:m)/m*2*pi); 0.01*cos((1:m)/m*2*pi); 0.01*sin((1:m)/m*2*pi); 0.1*sin((1:m)/m*2*pi)];
reference_vary = zeros(4, m);
for i=1:m
    reference_vary(:,i) = reference + varyRef(:,i);
end
simQuad( sys, innerController, x0, T, reference_vary,filter);

%% Final simulation
fprintf('Running the FINAL NL model simulation...\n')
close all
sim('simulation2.mdl')

%% BONUS - Lemniscate model
fprintf('BONUS - Lemniscate model...\n')
close all
sim('simulation_BONUS.mdl')

%% BONUS - Slew rate constraints
% run after doing nonlinear simulations otherwise the NL simulations won't
% work (because of the additional controller argument)
fprintf('BONUS - SLEW RATE CONSTRAINTS...\n')
yalmip('clear');
close all;
Q = diag([250, 1000, 1000, 5, 0.1, 0.1, 1]);% corresponding to first poles
% Q = diag([250, 10000, 10000, 5, 0.1, 0.1, 5]);
R = eye(4);
[K,P] = dlqr(sys.A,sys.B,Q,R); %K=-Klqr
% Declare variables
x = sdpvar(7,N,'full');
u = sdpvar(4,N,'full');
ref = sdpvar(4,1,'full');
steady_x = sdpvar(7,1,'full');
steady_u = sdpvar(4,1,'full');
dest = sdpvar(7,1,'full');
u_prev = sdpvar(4,N,'full'); 
tolerance = 0.003;

% Define constraints and objective
con = [];
con = con + (steady_mat*[steady_x;steady_u] ==  [Bd*dest; ref-Cd*dest]);
con = con + (x(:,2) == sys.A*x(:,1) + sys.B*u(:,1) + dest); % System dynamics
con = con + (-zDotMax<=x(1,:)<=zDotMax);
% con = con + (-alphaMax<=x(2,:)<=alphaMax); 
% con = con + (-betaMax<=x(3,:)<=betaMax);
con = con + (-alphaDotMax<=x(5,:)<=alphaDotMax);
con = con + (-betaDotMax<=x(6,:)<=betaDotMax);
con = con + (-gammaDotMax<=x(7,:)<=gammaDotMax);
con = con + (uMin<=u<=uMax);
con = con + (uMin<=u_prev<=uMax);
con = con + (abs(u(:,1)- u_prev(:,1))<= tolerance);

obj=(x(:,1)-steady_x)'*Q*(x(:,1)-steady_x) + (u(:,1)-steady_u)'*R*(u(:,1)-steady_u);
for i = 2:(N-1)
    con = con + (u(:,i-1) == u_prev(:,i));
    con = con + (abs(u(:,i)- u_prev(:,i))<= tolerance);  
    con = con + (x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i) + dest); % System dynamics
    obj = obj + (x(:,i)-steady_x)'*Q*(x(:,i)-steady_x) + (u(:,i)-steady_u)'*R*(u(:,i)-steady_u); % Cost function
end
obj = obj + (x(:,N)-steady_x)'*P*(x(:,N)-steady_x); % Terminal weight

% Compute optimizer
innerController = optimizer(con, obj, [], [x(:,1);ref;dest], u(:,1));

%Constant reference
reference = [0.4; 2*d2r; -2*d2r; 60*d2r];
simQuad(sys, innerController, x0, T, reference,filter);
%Varying reference
m=T/sys.Ts;
ref1 = [0.4; 2*d2r; -2*d2r; 60*d2r];
varyRef = [0.05*sin((1:m)/m*2*pi); 0.01*cos((1:m)/m*2*pi); 0.01*sin((1:m)/m*2*pi); 0.1*sin((1:m)/m*2*pi)];
reference_vary = zeros(4, m);
for i=1:m
    reference_vary(:,i) = ref1 + varyRef(:,i);
end
% simQuad( sys, innerController, x0, T, reference_vary,filter);






