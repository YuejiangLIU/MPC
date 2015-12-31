% MPC Ex5 | Yuejiang LIU | Apr. 2015
function ex5()
clear all; clc; close all;

PLOTOUT = 1;

nx = 2;
nd = 1;
A = [0.7115 -0.4345; 0.4345 0.8853];
B = [0.2173; 0.0573];
C = [0 1];
Bd = zeros(nx,1);
Cd = 1;

AA = [A Bd; zeros(nd,nx) eye(nd)];
CC = [C Cd];
p = [0.5 0.6 0.7];
L = place(AA',-CC',p).';
return;
% initial state
d = 0.2;
u0 = 0;
x(:,1) = [1;2];
est.x(:,1) = [3;0];
est.d(1) = 0;
T = 50;

% reference 
r = 0.5;

% MPC parameter
N = 5;
Q = eye(nx); 
R = 1;
P = dlyap(A,Q);

%% --- steady state controller ---
% Define optimization variables
xs = sdpvar(2,1,'full');
us = sdpvar(1,1,'full');
ds = sdpvar(1,1,'full');
% Define constraints and objective
con = [];
con = con + [[eye(nx)-A -B;C 0]*[xs;us] == [Bd*ds;r-Cd*ds]];
con = con + [-3<= us <= 3];
obj = us'*us;
% Compile the matrices
option = sdpsettings('solver','mosek');
%sol = optimize(con,obj,option);
optSteady = optimizer(con, obj, option, ds, [xs;us]);

%% --- tracking MPC controller ---
% Define optimization variables
dx = sdpvar(2,N,'full');
du = sdpvar(1,N,'full');
us = sdpvar(1,1,'full');
% Define constraints and objective
con = [];
obj = 0;
for i = 1:N-1
con = con + [dx(:,i+1) == A*dx(:,i) + B*du(:,i)]; % System dynamics
con = con + [-3 <= du(:,i)+us <= 3]; % Input constraints
obj = obj + dx(:,i)'*Q*dx(:,i) + du(:,i)'*R*du(:,i); % Cost function
end
obj = obj + dx(:,N)'*P*dx(:,N); % Terminal weight
% Compile the matrices
option = sdpsettings('solver','mosek');
ctrl = optimizer(con, obj, option, [dx(:,1);us], du(:,1)+us);

%% --- disturbance estimation ---
for tt = 1:T
    % steady state
    steady(:,tt) = optSteady{est.d(tt)};
    st.x(:,tt) = steady(1:2,tt);
    st.u(1,tt) = steady(3,tt);
    st.y(1,tt) = C*st.x(:,tt) + Cd*est.d(tt);
    % sys simulation
    x(:,tt+1) = syst(A,B,x(:,tt),u0);
    y(tt) = obs(C,d,x(:,tt));
    % estimate simulation
    [est.x(:,tt+1),est.d(tt+1)] = estimate(0,y(tt),A,B,C,Bd,Cd,L,est.x(:,tt),est.d(tt));
end

%% --- estimation visualization --- 
if ~PLOTOUT
    % plot state
    figure;
    t = 1:T+1;
    plot(t,x,t,est.x);
    legend('real x1','real x2','estimated x1','estimated x2');
    title('State Estimation');
    % plot disturbance
    figure;
    plot(t,est.d);
    ylim([-1 1]);
    title('Disturbance Estimation');
    % plot steady state 
    %figure;
    %t = 1:T;
    %plot(t,st.x,t,st.u,'-o',t,st.y,'-*');
    %ylim([-5 5]);
    %title('Steady State');  
    % plot output
    %figure;
    %plot(t,y);
    %title('Output');
end

%% --- tracking simulation ---
for tt = 1:T    
    % steady state target
    steady(:,tt) = optSteady{est.d(tt)};    % steady optimizer
    st.x(:,tt) = steady(1:2,tt);
    st.u(1,tt) = steady(3,tt);  
    
    % MPC tracking 
    u(1,tt) = ctrl{[est.x(:,tt)-st.x(:,tt);st.u(1,tt)]};
    
    % real observation 
    y(tt) = obs(C,d,x(:,tt));
    % next state estimation
    [est.x(:,tt+1),est.d(tt+1)] = estimate(u(1,tt),y(tt),A,B,C,Bd,Cd,L,est.x(:,tt),est.d(tt));
    % next state real dynamics
    x(:,tt+1) = syst(A,B,x(:,tt),u(1,tt));
    
end

figure;
plot(x(1,:),x(2,:),est.x(1,:),est.x(2,:));
xlabel('x_1');
ylabel('x_2');
legend('real','estimation');
title(['State (r = ' num2str(r) ') | Yuejiang LIU']);
grid on;

figure;
t = 1:T;
subplot(2,1,1);
plot(t,y);
xlabel('step');
ylabel('output');
title(['Tracking (r = ' num2str(r) ') | Yuejiang LIU']);
grid on;
subplot(2,1,2);
plot(t,u);
ylim([-4 4]);
xlabel('step');
ylabel('input');
grid on;

figure;
t = 1:T+1;
plot(t,d*ones(1,T+1));
hold on;
plot(t,est.d);
xlabel('step');
ylabel('disturbance');
legend('real','estimation');
title(['Disturbance (r = ' num2str(r) ') | Yuejiang LIU']);
grid on;

disp('Done');

end

function x = syst(A,B,x,u)
x = A*x + B*u;
end

function y = obs(C,d,x)
y = C*x + d;
end

function [x,d] = estimate(u,y,A,B,C,Bd,Cd,L,x,d)
[nx,nd] = size(B);
estimation = [ A Bd; zeros(nd,nx) eye(nd)] * [x;d] + [B;0] * u + L*(C*x + Cd*d - y);
x = estimation(1:nx);
d = estimation(1+nx:end);
end
