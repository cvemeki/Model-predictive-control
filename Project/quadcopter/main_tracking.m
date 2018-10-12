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

%% %%%%%%%%%%%%%% First MPC controller %%%%%%%%%%%%%%%%%%%
%initialization
d2r = 2*pi/360;%convert deg to rad
T = 10; %simulation time
N=2/sys.Ts; %horizon length

Q = diag([50, 10000, 10000, 2, 0.1, 0.1, 0.1]);
R = eye(4);
%definition of constraints
c1 = 1;
c2 = d2r*10;
c3 = d2r*10;
c5 = d2r*15;
c6 = d2r*15;
c7 = d2r*60;
%definition of the initial state Z_rate alpha beta gamma alpha_rate beta_rate gamma_rate 
x0 = [-1 10*d2r -10*d2r 120*d2r 0 0 0]';
% H*x<=h
Hx = [-1 0 0 0 0 0 0;1 0 0 0 0 0 0;
      0 -1 0 0 0 0 0;0 1 0 0 0 0 0;
      0 0 -1 0 0 0 0;0 0 1 0 0 0 0;
      0 0 0 0 0 0 0;0 0 0 0 0 0 0;%gamma no constraint
      0 0 0 0 -1 0 0;0 0 0 0 1 0 0;
      0 0 0 0 0 -1 0;0 0 0 0 0 1 0;
      0 0 0 0 0 0 -1;0 0 0 0 0 0 1;];
hx = [c1;c1;c2;c2;c3;c3;0;0;c5;c5;c6;c6;c7;c7;];
% Hu*u<=hu
Hu = [-1 0 0 0;
      1 0 0 0;
      0 -1 0 0;
      0 1 0 0;
      0 0 -1 0;
      0 0 1 0;
      0 0 0 -1;
      0 0 0 1;];
hu = [0;1;0;1;0;1;0;1;]-us(1)*[-1;1;-1;1;-1;1;-1;1;];
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
%initialization
T = 10; %simulation time
N=2/sys.Ts; %horizon length
Q = diag([50, 10000, 10000, 2, 0.1, 0.1, 0.1]);
R = eye(4);
% steady state computation
% Compute matrix for steady state
reference = [0.7; 9*d2r; 9*d2r; 40*d2r];
C = [eye(4) zeros(4,3)];
steady_mat = [eye(7)-sys.A -sys.B; C zeros(4,4)];
% result = (steady_mat)^-1*[zeros(7,1); reference];
% x_r = result(1:7)
% u_r = result(8:end)

% Define optimization variables
x = sdpvar(7,N,'full');
u = sdpvar(4,N,'full');
ref = sdpvar(4,1,'full');
steady = sdpvar(11,1,'full');

% Define constraints and objective
con = [];
con = con + (steady ==  (steady_mat)^-1*[zeros(7,1); ref]);
con = [con, x(:,2)-steady(1:7) == sys.A*(x(:,1)-steady(1:7)) + sys.B*(u(:,1)-steady(8:end))]; % System dynamics
con = con + (abs(x(1,:)) <= c1); % State constraints
con = con + (abs(x(2,:)) <= c2);
con = con + (abs(x(3,:)) <= c3);
con = con + (abs(x(5,:)) <= c5);
con = con + (abs(x(6,:)) <= c6);
con = con + (abs(x(7,:)) <= c7);
con = con + (Hu*u(:,1) <= hu); % Input constraints

obj=(x(:,1)-steady(1:7))'*Q*(x(:,1)-steady(1:7)) + (u(:,1)-steady(8:end))'*R*(u(:,1)-steady(8:end));
for i = 2:(N-1)
    con = [con, x(:,i+1)-steady(1:7) == sys.A*(x(:,i)-steady(1:7)) + sys.B*(u(:,i)-steady(8:end))]; % System dynamics
    con = [con, Hu*u(:,i) <= hu];
    obj = obj + (x(:,i)-steady(1:7))'*Q*(x(:,i)-steady(1:7)) + (u(:,i)-steady(8:end))'*R*(u(:,i)-steady(8:end)); % Cost function
end
obj = obj + (x(:,N)-steady(1:7))'*P*(x(:,N)-steady(1:7)); % Terminal weight

% Compute optimizer
innerController = optimizer(con, obj, [], [x(:,1);ref], u(:,1));
simQuad(sys, innerController, x0, T, reference);


%% Nonlinear model simulation - no disturbance
fprintf('Running the FIRST NL model simulation...\n')
close all
sim('simulation1.mdl') 


pause

%% Disturbance estimation
%estimator
A_est = [sys.A eye(length(sys.A));       %% A chapeau
         zeros(length(sys.A)) eye(length(sys.A))];
B_est = [sys.B; zeros(length(sys.A),4)]; %% B chapeau
C_est = [eye(length(sys.A)) zeros(length(sys.A))];
% L_x = diag([0.5,0.5,0.5,0.5,0.5,0.5,0.5]);
% L_d = diag([0.5,0.5,0.5,0.5,0.5,0.5,0.5]);
% L = [L_x; L_d];

pole = 0.8:0.01:0.93;
% pole(2)=0.23;
% pole(3)=0.26;
% pole(11) = 0.5;
L = place(A_est',C_est', pole)'; % Observer
filter.Af = A_est-L*C_est;
filter.Bf = [B_est L];
%% Offset free MPC
fprintf('PART III - OFFSET FREE / Disturbance rejection...\n')
close all
%initialization
T = 10; %simulation time
N=3/sys.Ts; %horizon length
Q = diag([200, 10000, 10000, 2, 45, 45, 1]);
R = eye(4);
sys.Bd = eye(7);
% steady state computation
% Compute matrix for steady state
reference = 0.2*[0.3; 2*d2r; 2*d2r; 40*d2r];
C = [eye(4) zeros(4,3)];
steady_mat_disturbance = [eye(7)-sys.A -sys.B; C zeros(4,4)];
% result = (steady_mat)^-1*[zeros(7,1); reference];
% x_r = result(1:7)
% u_r = result(8:end)

% Define optimization variables
x = sdpvar(7,N,'full');
u = sdpvar(4,N,'full');
ref = sdpvar(4,1,'full');
steady = sdpvar(11,1,'full');
dest = sdpvar(7,1,'full');

% Define constraints and objective
con = [];
con = con + (steady ==  (steady_mat)^-1*[sys.Bd*dest; ref - zeros(4,7)*dest]);
con = [con, x(:,2)-steady(1:7) == sys.A*(x(:,1)-steady(1:7)) + sys.B*(u(:,1)-steady(8:end))]; % System dynamics
% con = con + (abs(x(1,:)) <= c1); % State constraints
% con = con + (abs(x(2,:)) <= c2);
% con = con + (abs(x(3,:)) <= c3);
% con = con + (abs(x(5,:)) <= c5);
% con = con + (abs(x(6,:)) <= c6);
% con = con + (abs(x(7,:)) <= c7);
% con = con + (Hu*u(:,1) <= hu); % Input constraints
con = con + (-0.7007 <= u <= 0.2993); % Input constraints

obj=(x(:,1)-steady(1:7))'*Q*(x(:,1)-steady(1:7)) + (u(:,1)-steady(8:end))'*R*(u(:,1)-steady(8:end));
for i = 2:(N-1)
    con = [con, x(:,i+1)-steady(1:7) == sys.A*(x(:,i)-steady(1:7)) + sys.B*(u(:,i)-steady(8:end))]; % System dynamics
%     con = [con, Hu*u(:,i) <= hu];
    obj = obj + (x(:,i)-steady(1:7))'*Q*(x(:,i)-steady(1:7)) + (u(:,i)-steady(8:end))'*R*(u(:,i)-steady(8:end)); % Cost function
end
% obj = objj + (x(:,N)-steady(1:7))'*P*(x(:,N)-steady(1:7)); % Terminal weight
% Compute optimizer
innerController = optimizer(con, obj, [], [x(:,1);ref;dest], u(:,1));
ref1 = reference*2;
ref2 = reference/3;
reference_vary = [repmat(ref1,1,5/sys.Ts) repmat(ref2,1,10/sys.Ts + 1)];
simQuad(sys, innerController, x0, 15, reference_vary,filter);
   



%% Final simulation
fprintf('Running the FINAL NL model simulation...\n')
sim('simulation2.mdl') 
pause
%% BONUS - Slew rate constraints
% run after doing nonlinear simulations otherwise the NL simulations won't
% work (because of the additional controller argument)
fprintf('BONUS - SLEW RATE CONSTRAINTS...\n')






