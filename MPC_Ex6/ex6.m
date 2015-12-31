%% Solve a box-constrained optimal control problem using ADMM
clear all;
clc;
close all;
randn('state',0);
rand('state',0);

% Set to 1 for Part 1 of the exercise and to 2 for Part 2
system = 2;

%% >>>>>>>>>>>>>>> Defining problem's parameters >>>>>>>>>>>>>>>>>%%
% System's dimensions and control horizon
% n - states, m - controls
if system == 1
    A = [1.988 -0.998; 1 0]; 
    B = [0.125; 0];
    n = size(A,1); m = size(B,2); N = 30; rho = 5; % ADMM dual stepsize
    x_init = 5*randn(n,1);
    Q = eye(n); R = 15*eye(m);
    mat = blkdiag(Q,R);
    [~,S,~] = dlqr(A,B,Q,R,zeros(n,m));
    % constraints
    umax = 10; 
    umin = -umax;
    xmax = 9;
    xmin = -xmax;
elseif system == 2
    n = 20; m = 5; N = 20; rho = 10; x_init = 4*randn(n,1);
    A = randn(n);
    A = A/max(abs(eig(A)));
    B = randn(n,m);
    B = 1.1*B / max(abs(svd(B)));
    mat = randn(n+m); mat = mat*mat';
    mat(1:n,n+1:end) = zeros(n,m);
    mat(n+1:end,1:n) = zeros(m,n);
    mat = sparse(mat);
    Q = full(mat(1:n,1:n));
    R = full(mat(n+1:end,n+1:end));
    [~,S,~] = dlqr(A,B,Q,R,zeros(n,m));
    % constraints
    umax = 4;
    umin = -umax;
    xmax = 9;
    xmin = -xmax;
end

%%>>>>>>>>>>>>>>> Yalmip <<<<<<<<<<<<<<<<%%
X = sdpvar(n,N+1);
U = sdpvar(m,N);
con = []; 
for i = 1:N
    con = con + (X(:,i+1) == A*X(:,i) + B*U(:,i));
end
con = con + (X(:,1) == x_init);
for i = 1:N
    con = con + (X(:,i+1) <= xmax);
    con = con + (X(:,i+1) >= xmin);
    con = con + (U(:,i) <= umax);
    con = con + (U(:,i) >= umin);
end
obj = 0;
for i = 1:N
    obj = obj + 0.5 * [X(:,i); U(:,i)]'*mat*[X(:,i); U(:,i)];
end
obj = obj + 0.5 * X(:,N+1)'*S*X(:,N+1);

ops = sdpsettings('solver','sedumi','verbose',0,'sedumi.eps',1e-12);
results = solvesdp(con, obj, ops);
t_yal = results.solvertime;
disp(sprintf('Yalmip (using Sedumi) took %f seconds to solve',t_yal))

%%>>>>>>>>>>>>>>> ADMM <<<<<<<<<<<<<<<<%%
if system == 1
    time_inst = [5 100 200 300]; % iters @ which I save the state trajectory
    tic
    [X1,X2,z,bar_z,history] = admm(Q, R, S, A, B, n, m, N, rho, x_init, umax, umin, xmax, xmin, time_inst);
    toc
    
    %% State trajectories plotting
    H = [eye(n); -eye(n)];  h = [ones(n,1)*xmax; -ones(n,1)*xmin];
    % Plot the state constraints
    O = Polyhedron(H,h);
    figure(1); clf;
    g=O.plot('color', 'yellow');
    hold on;
    
    colors = {'g','r','b','k'};
    for iii = 1:size(X1,1)
        h(iii) = plot(X1(iii,:),X2(iii,:),colors{iii},'LineWidth', 1.5);
        legend(h(iii), ['k = ',num2str(time_inst(iii))]);
    end
    ylmpf = plot(double(X(1,:)),double(X(2,:)),'k--');
    leg = legend([g;h;ylmpf],'State Constraints','k=5','k=100','k=200','k=300','yalmip');
    title('State trajectories evolution as ADMM converges')
    
else
    time_inst = [5 100 200 300]; % iters @ which I save the state trajectory
    tic
    [~,~,z,bar_z,history] = admm(Q, R, S, A, B, n, m, N, rho, x_init, umax, umin, xmax, xmin, time_inst);
    toc
end


%% Reporting
K = length(history.objval);

figure;
plot(1:K, history.objval, 'k', [1 K], [double(obj) double(obj)], 'k--', 'MarkerSize', 10, 'LineWidth', 2);
l1 = ylabel('$V(x,u)$', 'fontsize', 14); 
xlabel('iter (k)');
legend('Objective','Optimal value');
set(l1,'Interpreter','LaTex');

% Comparison with Yalmip
x = z(1:n*(N+1));   u = z(n*(N+1)+1:end);
x = reshape(x,n,N+1); u = reshape(u,m,N);
fprintf('                OPTIMAL SOLUTION from Yalmip               \n')
fprintf('x = \n'); disp(double(X(:,1:10)))
fprintf('u = \n'); disp(double(U(:,1:10)))
fprintf('                OPTIMAL SOLUTION from ADMM               \n')
fprintf('x = \n'); disp(x(:,1:10))
fprintf('u = \n'); disp(u(:,1:10))

if system ==1
    figure;
    plot(1:N, u(:,1:N), 'k', 1:N, double(U(:,1:N)), 'k--', 'MarkerSize', 10, 'LineWidth', 2);
    l1 = ylabel('Input', 'fontsize', 14);
    xlabel('iter (i)');
    legend('ADMM','yalmip');
end

fprintf('                ERROR TO YALMIP for ADMM              \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x-double(X)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u-double(U)));