% MPC Ex3 Yuejiang LIU
clc;clear all;close all;

%% Q1
alpha = pi/6;
beta = 0.8;
A = [ cos(alpha) sin(alpha);-sin(alpha) cos(alpha)]*beta;
gamma = pi/3;
H = [ cos(gamma) sin(gamma);-cos(gamma) -sin(gamma);sin(gamma) -cos(gamma); -sin(gamma) cos(gamma)];
h = [2; 1; 2; 5];

Fnow = H;
fnow = h;
figure;
Pori = Polyhedron(Fnow,fnow);

while 1
    Pnow = Polyhedron(Fnow,fnow);
    Fpre = Fnow*A;
    fpre = fnow;
    preP = Polyhedron(Fpre,fpre);
    Fnext = [Fpre;Fnow];
    fnext = [fpre;fnow];
    Pnext = Polyhedron(Fnext,fnext);
    if Pnext==Pnow 
        Pmax = Pnext;
        break;
    else
        %plot([Pori Pnext]);
        %xlabel('x_1'); ylabel('x_2');
        %clf;
        Fnow = Fnext;
        fnow = fnext;
    end
end
plot([Pori Pmax]);
legend('State Constraint Set','Maximal Invariant Set');
xlabel('x_1'); ylabel('x_2');
title('Maximal Invariant Set (Yuejiang LIU)');

B = [0;0];
C = eye(2);
D = 0;
dt = 1;
t = 0:dt:30;
u = zeros(size(t));
sys = ss(A,B,C,D,dt);
x0 = [-3;3];
[y1 t x1] = lsim(sys,u,t,x0);

x0 = [-3;1];
[y2 t x2] = lsim(sys,u,t,x0);
x0 = [-2;1];
[y3 t x3] = lsim(sys,u,t,x0);
x0 = [-1;1];
[y4 t x4] = lsim(sys,u,t,x0);
x0 = [1;1];
[y5 t x5] = lsim(sys,u,t,x0);
x0 = [2;1];
[y6 t x6] = lsim(sys,u,t,x0);
hold on; 
plot(x1(:,1),x1(:,2),'b:o',...
    x2(:,1),x2(:,2),'y:*',...
    x3(:,1),x3(:,2),'y:s',...
    x4(:,1),x4(:,2),'y:+',...
    x5(:,1),x5(:,2),'y:x',...
    x6(:,1),x6(:,2),'y:d');

%% Q2
alpha = pi/6;
beta = 0.8;
A = [ cos(alpha) sin(alpha);-sin(alpha) cos(alpha)]*beta;
gamma = pi/3;
H = [ cos(gamma) sin(gamma);-cos(gamma) -sin(gamma);sin(gamma) -cos(gamma); -sin(gamma) cos(gamma)];
h = [2; 1; 2; 5];

B = [0.5;0.5];
G = [-1;1];
g = [0.5;0.5];

Fnow = H;
fnow = h;
Pori = Polyhedron(Fnow,fnow);

n = length(A);
q = length(G);

while 1
    Pnow = Polyhedron(Fnow,fnow);
    preP = projection( Polyhedron([Fnow*A Fnow*B;zeros(q,n) G], [fnow;g]) , [1:n]);
    Fpre = preP.A;
    fpre = preP.b;
    Fnext = [Fpre;Fnow];
    fnext = [fpre;fnow];
    Pnext = Polyhedron(Fnext,fnext);
    if Pnext==Pnow 
        Pmax = Pnext;
        break;
    else
        %plot([Pori Pnext]);
        %xlabel('x_1'); ylabel('x_2');
        %clf;
        Fnow = Fnext;
        fnow = fnext;
    end
end


Q = eye(2);
R = 1;
k = dlqr(A,B,Q,R);
AA = A-B*k;
Fnow = [H;G*k];
fnow = [h;g];

sys = ss(AA,B,C,D,dt);

while 1
    Pnow = Polyhedron(Fnow,fnow);
    Fpre = Fnow*AA;
    fpre = fnow;
    preP = Polyhedron(Fpre,fpre);
    Fnext = [Fpre;Fnow];
    fnext = [fpre;fnow];
    Pnext = Polyhedron(Fnext,fnext);
    if Pnext==Pnow
        Plqr = Pnext;
        break;
    else
        %plot([Pori Pnext]);
        %xlabel('x_1'); ylabel('x_2');
        %clf;
        Fnow = Fnext;
        fnow = fnext;
    end
end
figure;
plot([Pori Pmax Plqr]);
legend('State Constraint Set','Maximal Control Invariant Set','Maximal LQR Invariant Set');
xlabel('x_1'); ylabel('x_2');
title('Maximal Control & LQR Invariant Set (Yuejiang LIU)');
