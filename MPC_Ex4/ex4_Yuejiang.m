% MPC Ex4 | Yuejiang LIU | 2015 April

function ex4()

clear all
clc;
close all;

%% Problem Setup
A = [0.9752 1.4544
    -0.0327 0.9315];
B = [0.0248;0.0327];
S_X = [ -1 0
        1 0
        0 -1
        0 1];
s_X = [5;5;0.2;0.2];
Sori = Polyhedron(S_X,s_X);
S_U = [-1;1];
s_U = [1.75;1.75];
Q = eye(length(A))*10;
R = 1;
       
%% Terminal LQR  
[k,P] = dlqr(A,B,Q,R);

AA = A-B*k;
Fnext = [S_X];
fnext = [s_X];

while 1
    Fnow = Fnext;
    fnow = fnext;
    Snow = Polyhedron(Fnow,fnow);
    Fpre = Fnow*AA;
    fpre = fnow;
    Fnext = [Fpre;Fnow];
    fnext = [fpre;fnow]; 
    Snext = Polyhedron(Fnext,fnext);
    if Snext==Snow
        Slqr = Snow;
        Flqr = Fnow;
        flqr = fnow;
        break;
    end
end

%% Terminal Validation

sys = LTISystem('A',A,'B',B);
sys.x.max = [5;0.2];
sys.x.min = [-5;-0.2];
sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

kk = sys.LQRGain;
PP = sys.LQRPenalty.weight;
SSlqr = sys.LQRSet;

if kk + k == 0
    fprintf('The terminal controller is correct.\n');
else
    fprintf('Warning: The terminal controller is not correct.\n');
end

if PP == P
    fprintf('The terminal weight is correct.\n');
else
    fprintf('Warning: The terminal weight is not correct.\n');
end

if SSlqr == Slqr
    fprintf('The terminal set is correct.\n');
else
    fprintf('Warning: The terminal set is not correct.\n');
end


%% YMP controller 

F = S_X;
f = s_X;
M = S_U;
m = s_U;
Ff = Flqr;
ff = flqr;
Qf = P;
N = 10;
% Define optimization variables
x = sdpvar(2,N,'full');
u = sdpvar(1,N,'full');
% Define constraints and objective
con = [];
obj = 0;
for i = 1:N-1
con = con + [x(:,i+1) == A*x(:,i) + B*u(:,i)]; % System dynamics
con = con + [F*x(:,i) <= f]; % State constraints
con = con + [M*u(:,i) <= m]; % Input constraints
obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); % Cost function
end
con = con + [Ff*x(:,N) <= ff]; % Terminal constraint
obj = obj + x(:,N)'*Qf*x(:,N); % Terminal weight
% Compile the matrices
option = sdpsettings('solver','sedumi','sedumi.eps',1e-12);
ctrl = optimizer(con, obj, option, x(:,1), u(:,1));


%% Simulation 
format long 
x0 = [3;0];
%x0 = [2.886925397404798;-0.149094334873512];
% u0 = -1.559459782064587;
xnow = x0;
N = 10;
T = 50;
for tt = 1:T
    [u0,infeasible] = mpc(xnow,N,A,B,S_X,s_X,S_U,s_U,Q,R,Flqr,flqr,P);    
    %[u0,infeasible] = ctrl{xnow};
    if infeasible == 1 
        fprintf('Problem is infeasible at step %d.\n',tt);
        error('Unsolved.')
        u0 = 0;
    end
    u(tt) = u0;
    xnow = A*xnow + B*u0;
    x(:,tt) = xnow; 
end

figure;
t = 0:tt;
plot(t,[x0(1) x(1,:)],'--o',t,[x0(2) x(2,:)],'--v',t,[0 u],'--s');
xlabel('Time (s)');
legend('Position','Velocity','Input');
ylim([-5 5]);
grid on;

figure;
plot([Sori Slqr]);
hold on;
plot([x0(1) x(1,:)],[x0(2) x(2,:)],'--o');
legend('State Constraint Set','LQR Invariant Set','Trajectory');
xlabel('x_1'); ylabel('x_2');
grid on;


function [u0,infeasible] = mpc(xnow,N,A,B,F,f,M,m,Q,R,Ff,ff,Qf)
% MPC controller 

dimX = length(A);

T1_shape = tril(triu(ones(N,N),-1),-1);
T1_mat = -A;
T1 = kron(T1_shape,T1_mat);
T2_shape = eye(N);
T2_mat = eye(dimX);
T2 = kron(T2_shape,T2_mat);
T3_shape = eye(N);
T3_mat = -B;
T3 = kron(T3_shape,T3_mat);
T_shape = [(T1_shape+T2_shape) T3_shape];
T = [(T1+T2) T3];

t_shape = zeros(N,1);
t_shape(1,1) = 1;
t_mat = A*xnow;
t = kron(t_shape,t_mat);

G = blkdiag(kron(eye(N-1),F),Ff,kron(eye(N),M)); 

g1_shape = ones(N-1,1);
g1_mat = f;
g1 = kron(g1_shape,g1_mat);
g2_shape = ones(1);
g2_mat = ff;
g2 = kron(g2_shape,g2_mat);
g3_shape = ones(N,1);
g3_mat = m;
g3 = kron(g3_shape,g3_mat);
g_shape = [g1_shape;g2_shape;g3_shape];
g = [g1;g2;g3];

H = blkdiag(kron(eye(N-1),Q),Qf,kron(eye(N),R));
h = zeros(length(H),1);

[zopt, fval, flag] = quadprog(H,h,G,g,T,t);
infeasible = (flag ~= 1);

try
    u0 = zopt(N*dimX + 1);
catch 
    disp('Infeasible');
    u0 = 0;
end
