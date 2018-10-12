
% loads:
%    hovering equilibrium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outer controller optimizer instance
load('quadData.mat')
outerController = getOuterController(Ac, 'cplex');%cplex
disp('Data successfully loaded')



%% %%%%%%%%%%%%%% First MPC controller %%%%%%%%%%%%%%%%%%%
N = 5;
x1 = sdpvar(7,N);  %% create a 7xN matrix
constraints = [];
angleMax = pi/3;
for i = 2:N
    constraints = constraints + (x(:,i+1) == A*x(:,i) + B*u(:,i));  %System dynamics
end 
R = Ac;
Q = Bc;
objective = sum(x'*R*x,2) + sum(u'*Q*x,2);
options = sdpsettings('solver','qpip');   %% what is qpip??
x0 = [-1,10,-10,120,0,0,0]';
innnerController = optimizer(constraint, objective, options, x(:,1), u(:,1));
simQuad(sys,innerController, x0, T);
pause

%% Reference tracking - no disturbance, no invariant sets
fprintf('PART II - reference tracking...\n')

pause

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






