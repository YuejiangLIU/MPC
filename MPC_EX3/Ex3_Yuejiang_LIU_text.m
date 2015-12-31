% MPC Ex3 Yuejiang LIU
clc;clear;close;

A = [ 1 1; 0 1];
B = [ 1; 0.5];
Q = eye(2);
R = 90;

k = dlqr(A,B,Q,R);
C = [1 0];
D = 0;
sys_cl = ss(A-B*k,B,C,D);
AA = A-B*k;
%
H = [ -1 0
       1 0
       0 -1
       0 1
       -k(1) -k(2)
       k(1) k(2)];
h = [5 5 10 10 0.1 0.1]';

%%
Fnow = H;
fnow = h;
figure;
Pori = Polyhedron(Fnow,fnow);

while 1
    Pnow = Polyhedron(Fnow,fnow);
    Fpre = Fnow*AA;
    fpre = fnow;
    preP = Polyhedron(Fpre,fpre);
    Fnext = [Fpre;Fnow];
    fnext = [fpre;fnow];
    Pnext = Polyhedron(Fnext,fnext);
    if Pnext==Pnow 
        break;
    else
        %plot([Pori Pnext]);
        %xlabel('x_1'); ylabel('x_2');
        %clf;
        Fnow = Fnext;
        fnow = fnext;
    end
end
plot([Pori Pnext]);
xlabel('x_1'); ylabel('x_2');
title('Maximum Invariant Set');
Fnow;
fnow;

