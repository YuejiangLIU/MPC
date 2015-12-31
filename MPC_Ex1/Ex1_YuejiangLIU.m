% MPC Ex1 | Yuejiang LIU
clc;
close all;
clear all;

A = [ 4/3 -2/3; 1 0 ];
B = [ 1; 0 ];
C = [ -2/3 1 ];
D = 0;


Q = C'*C + 0.001*eye(2);
R =  0.001;
Pf = Q;

x0 = [10;10];

dt = 1;
t = 0:dt:30;
u = zeros(size(t));


sys = ss(A,B,C,D,dt);

%{
sysc = d2c(sys);

figure;
impulse(sys);
grid on;
figure;
step(sys,sysc);
grid on;

figure;
lsim(sys,u,t,x0);
grid on;
%}

%% Ex1


%% Ex2

N = 6;
dN = 1;
for i = 1:4
    [K0,S0] = finite_dlqr(A,B,Q,R,Pf,N);
    sys_cl = ss(A+B*K0,B,C,D,dt);
    [y(:,i) t x(:,:,i)] = lsim(sys_cl,u,t,x0);
    M{i} = A'*S0*A - S0;
    legendinfo{i} = ['N = ' num2str(N)];
    N = N + dN;
end

figure;
for i = 1:4
    plot(x(:,1,i),x(:,2,i),'-o');
    axis([-20 20 -20 20]);
    hold on;
end
hold off;
xlabel('x_1');
ylabel('x_2');
grid on;
legend(legendinfo);
title('Closed-loop Response (a)');

figure;
plot(t,y(:,1:end));
ylim([-2 10]);
hold on;




%% Ex3

[K,S] = dlqr(A,B,Q,R);
sys_cl = ss(A-B*K,B,C,D,dt);

[y t] = lsim(sys_cl,u,t,x0);
plot(t,y,'k');
legendinfo{i+1} = 'LQR';
legend(legendinfo);
grid on;
title('Closed-loop Response (b)');
xlabel('Time (s)');
ylabel('Output');


format long
cost_infinite = x0'*S*x0;
[K0,S0] = finite_dlqr(A,B,Q,R,Pf,8);
cost_finite = x0'*S0*x0;
[cost_infinite cost_finite cost_infinite-cost_finite] 
% test = lqr(sysd,Q,R)
A'*S0*A - S0;
A'*S*A - S;
