clear all;
close all;
clc;
% loads:
%    hovering equilibrium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outer controller optimizer instance
load('quadData.mat')
%outerController = getOuterController(Ac, 'cplex');
outerController = getOuterController(Ac,'gurobi');
disp('Data successfully loaded')
pause;
%% %%%%%%%%%%%%%% First MPC controller %%%%%%%%%%%%%%%%%%%
close all;
fprintf('PART I - MPC controller...\n')
% Parameter 
Lx = size(sys.A,2);
Lu = size(sys.B,2);
Q = diag([50 1e4 1e4 2 1e-1 1e-1 1e-1]);
%Q = diag([50 1e4 1e4 2 1 1 10]);
%Q = diag([5 1e3 1e3 5 1e-1 1e-1 1e-1]); 
R = eye(Lu)*1;
%[k,Qf] = dlqr(sys.A,sys.B,Q,R);
Npred = 20;
% Constraint Constant 
zdotMax = 1;
angleMax = 10/180*pi;
angledotMax = 15/180*pi;
yawdotMax = 60/180*pi;
uMin  = -us;
uMax  = ones(4,1)-us;
yawMax = 5;
% Terminal Constraint
syst = LTISystem('A',sys.A,'B',sys.B);
syst.x.max = [zdotMax;angleMax;angleMax;yawMax;angledotMax;angledotMax;yawdotMax];
syst.x.min = [-zdotMax;-angleMax;-angleMax;-yawMax;-angledotMax;-angledotMax;-yawdotMax];
syst.x.penalty = QuadFunction(Q);
syst.u.penalty = QuadFunction(R);
Qf = syst.LQRPenalty.weight;
Sf = syst.LQRSet;
Ff = Sf.A;
ff = Sf.b;

% --- Controller Start ---  
% Define optimization variables
x = sdpvar(Lx,Npred,'full');
u = sdpvar(Lu,Npred,'full');
% Define constraints and objective
con = [];
obj = 0;
for i = 1:Npred-1
con = con + [x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i)]; % System dynamics
con = con + [-zdotMax <= x(1,i) <= zdotMax]; % State constraints
con = con + [-angleMax <= x(2:3,i) <= angleMax];
con = con + [-angledotMax <= x(5:6,i) <= angledotMax];
con = con + [-yawdotMax <= x(7,i) <= yawdotMax];
con = con + [uMin <= u(:,i) <= uMax]; % Input constraints
obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); % Cost function
end
con = con + [Ff * x(:,Npred) <= ff]; % Terminal constraint
obj = obj + x(:,Npred)'*Qf*x(:,Npred); % Terminal weight
% Compile the matrices
option = sdpsettings('solver','gurobi');
innerController = optimizer(con, obj, option, x(:,1), u(:,1));
% --- Controller End ---  

zdot0 = -1;
roll0 = 10/180*pi;
pitch0 = -10/180*pi;
yaw0 = 120/180*pi;
x0 = [zdot0;roll0;pitch0;yaw0;0;0;0]*0.9;
T = 3;
simQuad(sys, innerController, x0, T);
pause;
%% Reference tracking - no disturbance, no invariant sets
fprintf('PART II - reference tracking...\n');
yalmip('clear');
close all;

%Q = diag([50 1e3 1e3 10 0 0 0]);
Q = diag([50 1e3 1e3 2 1e-1 1e-1 1e-1]);
R = eye(Lu)*1;
Npred = 20;

C = [eye(4) zeros(4,3)];
% vvvvvvvvvvvvvvv steady state vvvvvvvvvvvvvvvvvvvvv
% Define optimization variables
xr = sdpvar(Lx,1,'full');
ur = sdpvar(Lu,1,'full');
ref = sdpvar(4,1,'full');
% Define constraints and objective
con = [];
con = con + [sys.A*xr + sys.B*ur == xr];
con = con + [C*xr == ref];
con = con + [-zdotMax <= xr(1) <= zdotMax]; % State constraints
con = con + [-angleMax <= xr(2:3) <= angleMax];
con = con + [-angledotMax <= xr(5:6) <= angledotMax];
con = con + [-yawdotMax <= xr(7) <= yawdotMax];
con = con + [uMin <= ur <= uMax]; % Input constraints
obj = ur'*ur;
% Compile the matrices
option = sdpsettings('solver','gurobi');
%sol = optimize(con,obj,option);
optSteady = optimizer(con, obj, option, ref, [xr;ur]);
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

zdotr = -1;
rollr = 10/180*pi;
pitchr = -10/180*pi;
yawr = 120/180*pi;
ref = [zdotr rollr pitchr yawr]'/2;
%output = optSteady{ref};                        % == [ref;us];
%xr = output(1:Lx,1)
%ur = output(Lx+1:end,1)

% vvvvvvv tracking MPC controller vvvvvvvvvv
% Define optimization variables
xc = sdpvar(Lx,1);
rc = sdpvar(4,1);
ur = zeros(Lu,1);
xr = [rc;0;0;0];
dx = sdpvar(Lx,Npred,'full');
du = sdpvar(Lu,Npred,'full');
% Define constraints and objective
con = [];
obj = 0;
for i = 1:Npred-1
con = con + [dx(:,1) == xc - xr];
con = con + [dx(:,i+1) == sys.A*dx(:,i) + sys.B*du(:,i)]; % System dynamics
con = con + [-zdotMax <= dx(1,i) + xr(1,1) <= zdotMax]; % State constraints
con = con + [-angleMax <= dx(2:3,i) + xr(2:3,1) <= angleMax];
con = con + [-angledotMax <= dx(5:6,i) + xr(5:6,1) <= angledotMax];
con = con + [-yawdotMax <= dx(7,i) + xr(7,1) <= yawdotMax];
con = con + [uMin <= du(:,i) <= uMax]; % Input constraints
obj = obj + dx(:,i)'*Q*dx(:,i) + du(:,i)'*R*du(:,i); % Cost function
end
%con = con + [Ff * (dx(:,Npred)+xr) <= ff]; % Terminal constraint
%con = con + [Ff * dx(:,Npred) <= ff]; % Terminal constraint
con = con + [dx(:,Npred) == 0];
obj = obj + dx(:,Npred)'*Qf*dx(:,Npred); % Terminal weight
% Compile the matrices
option = sdpsettings('solver','gurobi');
innerController = optimizer(con, obj, option, [xc;rc], du(:,1));
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% --- constant reference tracking ---
zdotr = -0.8;
rollr = 8/180*pi;
pitchr = -6/180*pi;
yawr = 3/180*pi;
ref = [zdotr rollr pitchr yawr]'/2;
x0 = zeros(Lx,1);
simQuad( sys, innerController, x0, T, ref);

% --- varying reference tracking ---
ref1 = ref*2;
ref2 = ref/3;
T = 15;
ref = [repmat(ref1,1,7/sys.Ts) repmat(ref2,1,8/sys.Ts + 1)];
simQuad( sys, innerController, x0, T, ref);

pause;
%% Nonlinear model simulation - no disturbance
close all;
fprintf('Running the FIRST NL model simulation...\n')
sim('simulation1.mdl') 



%% Disturbance estimation
% Estimator
A_hat = [sys.A eye(size(sys.A,1))
         zeros(size(sys.A)) eye(size(sys.A,1))];
B_hat = [sys.B 
        zeros(size(sys.A,1),size(sys.B,2))];
C_hat = [eye(7) zeros(7)];

pole = [0.3 % d_yawdot
        0.5 % x_yawdot
        0.5 % x_pitchdot
        0.4 % d_yaw
        0.5 % x_yaw
        0.45 % d_roll
        0.45 % d_ptich 
        0.45 % d_pitchdot 
        0.4 % d_zdot 
        0.5 % x_rolldot 
        0.4 % d_rolldot 
        0.5 % x_pitch 
        0.5 % x_roll
        0.5 % x_zdot
        ];
        
Q = diag([10 1e3 1e3 5 1e-1 1e-1 1e-1]); 
%pole = [0.8 ones(1,6) zeros(1,7)]'
%pole = 0.4:0.02:0.02*13+0.4;
%poleend = 0.8;
%pole = poleend:-0.03:-0.03*13+poleend;
%pole = [0.54;0.55;0.56;0.71;0.72;0.73;0.74;...
%        0.54;0.55;0.56;0.71;0.72;0.73;0.74];        % [dyaw, dpitch, drow, dz, ... ] 
L = place(A_hat',C_hat',pole)';
filter.Af = A_hat - L*C_hat;
filter.Bf = [B_hat L];

% Offset free MPC
close all;
fprintf('PART III - OFFSET FREE / Disturbance rejection...\n')
Npred = 20;
% vvvvvvv offest free tracking MPC controller vvvvvvvvvv
% Define optimization variables
xc = sdpvar(Lx,1);
rc = sdpvar(4,1);
ur = zeros(Lu,1);
xr = [rc;0;0;0];
d_est = sdpvar(Lx,1);
dx = sdpvar(Lx,Npred,'full');
du = sdpvar(Lu,Npred,'full');
% Define constraints and objective
con = [];
obj = 0;
for i = 1:Npred-1
con = con + [dx(:,1) == xc - xr];
%con = con + [[dx(:,i+1)+xr(:,1);d_est] == filter.Af*[dx(:,i)+xr(:,1);d_est] + filter.Bf*[du(:,i);dx(:,i+1)+xr(:,1)]]; % System dynamics xf = Af*xf + Bf*[u; x];
con = con + [dx(:,i+1) == sys.A*dx(:,i) + sys.B*du(:,i) + d_est];
con = con + [-zdotMax <= dx(1,i) + xr(1,1) <= zdotMax]; % State constraints
con = con + [-angleMax <= dx(2:3,i) + xr(2:3,1) <= angleMax];
con = con + [-angledotMax <= dx(5:6,i) + xr(5:6,1) <= angledotMax];
con = con + [-yawdotMax <= dx(7,i) + xr(7,1) <= yawdotMax];
con = con + [uMin <= du(:,i) <= uMax]; % Input constraints
obj = obj + dx(:,i)'*Q*dx(:,i) + du(:,i)'*R*du(:,i); % Cost function
end
%con = con + [dx(:,Npred) == 0]; % Terminal constraint
obj = obj + dx(:,Npred)'*Qf*dx(:,Npred); % Terminal weight
% Compile the matrices
option = sdpsettings('solver','gurobi');
innerController = optimizer(con, obj, option, [xc;rc;d_est], du(:,1));
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


ref = [zdotr rollr pitchr yawr]'/2;
simQuad( sys, innerController, x0, T, ref, filter);

ref1 = ref*1.5;
ref2 = ref/3;
T = 15;
ref = [repmat(ref1,1,7/sys.Ts) repmat(ref2,1,8/sys.Ts + 1)];
simQuad( sys, innerController, x0, T, ref, filter);

pause;

%% Final simulation
close all;
fprintf('Running the FINAL NL model simulation...\n')
sim('simulation2.mdl') 

%% BONUS - Slew rate constraints
% run after doing nonlinear simulations otherwise the NL simulations won't
% work (because of the additional controller argument)
fprintf('BONUS - SLEW RATE CONSTRAINTS...\n')
yalmip('clear');
close all;

slewMax = 5e-2;

% Estimator
A_hat = [sys.A eye(size(sys.A,1))
         zeros(size(sys.A)) eye(size(sys.A,1))];
B_hat = [sys.B 
        zeros(size(sys.A,1),size(sys.B,2))];
C_hat = [eye(7) zeros(7)];

pole = [0.3 % d_yawdot
        0.5 % x_yawdot
        0.5 % x_pitchdot
        0.4 % d_yaw
        0.5 % x_yaw
        0.45 % d_roll
        0.45 % d_ptich 
        0.45 % d_pitchdot 
        0.4 % d_zdot 
        0.5 % x_rolldot 
        0.4 % d_rolldot 
        0.5 % x_pitch 
        0.5 % x_roll
        0.5 % x_zdot
        ];
        
Q = diag([10 1e3 1e3 5 1e-1 1e-1 1e-1]); 

L = place(A_hat',C_hat',pole)';
filter.Af = A_hat - L*C_hat;
filter.Bf = [B_hat L];

Npred = 20;
% vvvvvvv offest free tracking MPC controller vvvvvvvvvv
% Define optimization variables
xc = sdpvar(Lx,1);
rc = sdpvar(4,1);
ur = zeros(Lu,1);
xr = [rc;0;0;0];
d_est = sdpvar(Lx,1);
dx = sdpvar(Lx,Npred,'full');
du = sdpvar(Lu,Npred,'full');
u_prev = sdpvar(Lu,1,'full');
% Define constraints and objective
con = [];
obj = 0;
con = con + [-slewMax <= du(:,1)-u_prev <= slewMax];
for i = 1:Npred-1
con = con + [dx(:,1) == xc - xr];
%con = con + [[dx(:,i+1)+xr(:,1);d_est] == filter.Af*[dx(:,i)+xr(:,1);d_est] + filter.Bf*[du(:,i);dx(:,i+1)+xr(:,1)]]; % System dynamics xf = Af*xf + Bf*[u; x];
con = con + [dx(:,i+1) == sys.A*dx(:,i) + sys.B*du(:,i) + d_est];
con = con + [-zdotMax <= dx(1,i) + xr(1,1) <= zdotMax]; % State constraints
con = con + [-angleMax <= dx(2:3,i) + xr(2:3,1) <= angleMax];
con = con + [-angledotMax <= dx(5:6,i) + xr(5:6,1) <= angledotMax];
con = con + [-yawdotMax <= dx(7,i) + xr(7,1) <= yawdotMax];
con = con + [uMin <= du(:,i) <= uMax]; % Input constraints
con = con + [-slewMax <= du(:,i+1)-du(:,i) <= slewMax]; 
obj = obj + dx(:,i)'*Q*dx(:,i) + du(:,i)'*R*du(:,i); % Cost function
end
%con = con + [dx(:,Npred) == 0]; % Terminal constraint
obj = obj + dx(:,Npred)'*Qf*dx(:,Npred); % Terminal weight
% Compile the matrices
option = sdpsettings('solver','gurobi');
innerController = optimizer(con, obj, option, [xc;rc;u_prev;d_est], du(:,1));
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% --- constant reference tracking ---
zdotr = -0.8;
rollr = 8/180*pi;
pitchr = -6/180*pi;
yawr = 3/180*pi;
ref = [zdotr rollr pitchr yawr]'/2;
x0 = zeros(Lx,1);

% --- varying reference tracking 2 ---
ref1 = zeros(4,1);
ref2 = ref;
T = 15;
ref = [repmat(ref1,1,9.5/sys.Ts) repmat(ref2,1,5.5/sys.Ts + 1)];
simQuad( sys, innerController, x0, T, ref, filter, [], 1);

