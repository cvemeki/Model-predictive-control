close all;


% loads:
%    hovering equilibr  ium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outer controller optimizer instance
load('quadData.mat')
outerController = getOuterController(Ac, 'cplex');
disp('Data successfully loaded')



%% %%%%%%%%%%%%%% First MPC controller %%%%%%%%%%%%%%%%%%%
%initialization
N=20; %horizon length
d2r = 2*pi/360;%convert deg to rad jk
T = 10; %simulation time
Q = 1*diag([50;10000;10000;2;0.1;0.1;0.1]); %state matrix
% R=[1 2 3 4;
%    2 3 4 1;
%    3 4 1 2;
%    4 1 2 3;]; %input matrix
% R_ = 100*rand(4);  
% R = tril(R,-1)+triu(R',0);  
R = 1*eye(4);
% R = 1*diag([10;0.01;0.01;0.8]);
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
hu = [0.70071429;1-0.70071429;0.70071429;1-0.70071429;0.70071429;1-0.70071429;0.70071429;1-0.70071429;];

sys.uMin  = -u(1)*ones(4,1);
sys.uMax  = u(2)*ones(4,1);
% sys.uMin  = -0*ones(4,1);
% sys.uMax  = ones(4,1);
sys.angleMin = -[1;1]*hx(2,1);
sys.angleMax =  [1;1]*hx(2,1);
sys.zVelMin = -hx(1,1);
sys.zVelMax = hx(1,1);
sys.angVelMin   = -hx(5:7,1);
sys.angVelMax   = hx(5:7,1);

%Define optimization variables  
x = sdpvar(7,N,'full');
u = sdpvar(4,N,'full');
%Define constraints and objective 
con = ( x(:,2) == sys.A*x(:,1) + sys.B*u(:,1)) + (Hu*u(:,1) <= hu); % constraint when n == 1
obj = x(:,1)'*Q*x(:,1)+u(:,1)'*R*u(:,1);
for i = 2:N-1
    con = con + (x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i));   %System dynamics
    con = con + (Hx*x(:,i) <= hx);                   %State constraints
    con = con + (Hu*u(:,i) <= hu);                   %Input constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); %Cost function
end
% con = con + (Xf.A*x(:,N) <= Xf.b);   %Terminal constraint
% con = con + (x(:,N) == 0);   %Terminal constraint % infeasible
% [K,P] = dlqr(sys.A,sys.B,Q,R); %K=-Klqr

% Xf = Terminalset(sys.A,sys.B,K,Hx,hx,Hu,hu);

% con = con + (Xf.A*x(:,N) <= Xf.b);   %Terminal constraint
% obj = obj + x(:,N)'*P*x(:,N);      %Terminal weight 

%Compile the matrices
options = sdpsettings('solver','quadprog');
innerController = optimizer(con,obj,options,x(:,1),u(:,1));
%Can now compute the optimal control input using
[xt, ut, t] = simQuad( sys, innerController, x0, T);


%% Reference tracking - no disturbance, no invariant sets
fprintf('PART II - reference tracking...\n')

ref = [0.4;0.1;0.1;0.9];
sys.C = [eye(4) zeros(4,3) ];
M_ss = [ eye(length(sys.A)) - sys.A, -sys.B; sys.C, zeros( size(sys.C,1),size(sys.B,2) ) ];
xu_ss = M_ss\[zeros(length(sys.A),1);ref]
xss = xu_ss(1:7,1);
uss = xu_ss(8:11,1);

d_x0 = x0 - xss;
% xu_ss_ = inv(M_ss)*[zeros(length(sys.A),1);ref]
hx_ss = hx - Hx*xss
hx
hu_ss = hu - Hu*uss

sys.uMin  = -hu_ss(1)*ones(4,1);
sys.uMax  = hu_ss(2)*ones(4,1);
% sys.uMin  = -0*ones(4,1);
% sys.uMax  = ones(4,1);
sys.angleMin = -[1;1]*hx_ss(2,1);
sys.angleMax =  [1;1]*hx_ss(2,1);
sys.zVelMin = -hx_ss(1,1);
sys.zVelMax = hx_ss(1,1);
sys.angVelMin   = -hx_ss(5:7,1);
sys.angVelMax   = hx_ss(5:7,1);

d_x = sdpvar(7,N,'full');
d_u = sdpvar(4,N,'full');
%Define constraints and objective 
con_ss = ( d_x(:,2) == sys.A*d_x(:,1) + sys.B*d_u(:,1)) + (Hu*d_u(:,1) <= hu_ss); % constraint when n == 1
obj_ss = d_x(:,1)'*Q*d_x(:,1)+d_u(:,1)'*R*d_u(:,1);
for i = 2:N-1
    con_ss = con_ss + (d_x(:,i+1) == sys.A*d_x(:,i) + sys.B*d_u(:,i));   %System dynamics
    con_ss = con_ss + (Hx*d_x(:,i) <= hx_ss);                   %State constraints
    con_ss = con_ss + (Hu*d_u(:,i) <= hu_ss);                   %Input constraints
    obj_ss = obj_ss + d_x(:,i)'*Q*d_x(:,i) + d_u(:,i)'*R*d_u(:,i); %Cost function
end

%Compile the matrices
options = sdpsettings('solver','quadprog');
innerController_ss = optimizer(con_ss,obj_ss,options,d_x(:,1),d_u(:,1));
%Can now compute the optimal control input using
[xt, ut, t] = simQuad(sys, innerController_ss, d_x0, T);
% pause 

x_real(:,1) = x0;
for i = 1 : 100
    u_real(:,i) = ut(:,i) + uss
    x_real(:,i+1) = sys.A*x_real(:,i) + sys.B*u_real(:,i);
end
y_real = C*x_real;
figure()
hold on
plot(t(1:100),u_real(1,:));
plot(t(1:101),y_real(1,:));
plot(t(1:101),y_real(2,:));
plot(t(1:101),y_real(3,:));
plot(t(1:101),y_real(4,:));
legend('u1', 'zdot','alpha','beta','gamma');
hold off
%% Nonlinear model simulation - no disturbance
fprintf('Running the FIRST NL model simulation...\n')
sim('simulation1.mdl') 

pause

%% Disturbance estimation
%estimator


%% Offset free MPC
fprintf('PART III - OFFSET FREE / Disturbance rejection...\n')

pause

%% Final simulation
fprintf('Running the FINAL NL model simulation...\n')
sim('simulation2.mdl') 
pause
%% BONUS - Slew rate constraints
% run after doing nonlinear simulations otherwise the NL simulations won't
% work (because of the additional controller argument)
fprintf('BONUS - SLEW RATE CONSTRAINTS...\n')


function terminal = Terminalset(A,B,K,Hx,hx,Hu,hu)
    F = A + B*K;
    Omega = Polyhedron([Hx;Hu*K],[hx;hu]);  %input constraints 
    i = 1;
    while (1)
        Omega_pre = Polyhedron(Omega.A*F,Omega.b);
        Omega_inter =  Polyhedron([Omega.A;Omega_pre.A],[Omega.b;Omega_pre.b]);
            if Omega_inter == Omega
                Omega_inf = Omega_inter;
                break;
            end
        if i <= 15
            Omega_inf = Omega_inter;
            break
        end
        i = i + 1
        Omega = Omega_inter;
    end
    terminal = Omega_inf;
end



