function [X1,X2,z,bar_z,history] = admm(Q, R, S, A, B, n, m, N, rho, x_init, umax, umin, xmax, xmin, time_inst)

%% run admm
    Mat = blkdiag(kron(eye(N), Q),S,kron(eye(N), R));
    QUIET    = 0;
    MAX_ITER = 3000;
    EPS_ABS   = 1e-4;
    EPS_REL   = 1e-4;  
    % Dimensions
    z = zeros(n*(N+1)+m*N,1);  % First primal variable - states and inputs
    bar_z = z; % Second primal variable - states and inputs
    mu = z;  % Dual variables
    sol = zeros(n*(N+1)+m*N,1);  % Vector of KKT solution

    if ~QUIET
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
          'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
    end
    
    %%*******************************************************************%%
    %% Build the linear system                                           %%
    %% You are given the matrices QQ, AA, b as presented in the slides   %% 
    %% Your QP reads                                                     %%  
    %% minimize (1/2)z'*QQ*z                                             %%
    %% subject to AAz = b                                                %%
    %%            l <= z <= u                                            %%
    %% Formulate the linear system that you will invert in the first step%% 
    %% of the algorithm.                                                 %% 
    %%*******************************************************************%%
    
    QQ = blkdiag(kron(eye(N),Q),S,kron(eye(N),R));
    QQ = QQ + rho*eye((N+1)*n+N*m);
    AA = kron(eye(N+1),eye(n)) + [zeros(n,n*(N+1));[kron(eye(N),-A)...
         zeros((N)*n,n)]];
    AA = [AA [zeros(n,m*N);kron(eye(N),-B)]];
    b = [x_init; zeros(n*N,1)];
    

    %%*************************** Your code here ************************%%
    %%                                                                   %%
    l = [ones(n*(N+1),1)*xmin;ones(m*N,1)*umin];
    u = [ones(n*(N+1),1)*xmax;ones(m*N,1)*umax];
    offlineMM = [QQ AA';AA zeros(size(AA,1),size(AA,1))];
    %%*******************************************************************%%
    
    % Activate for keeping the history of the state sequence 
    plot_traj = 1;
    X1 = [];  X2 = [];
    iii = 1;
    
    for k=1:MAX_ITER
        
       %% z - update : Solve the linear system                           %%
       %% You need to recover only z from the solution of the linear     %%
       %% system solve, where z = (x,u).                                 %%

       %%************************ Your code here ************************%%
       %%                                                                %%
       zz = offlineMM\[rho*(bar_z-mu);b];
       z = zz(1:length(z));
       %%****************************************************************%%
       
       if (plot_traj == 1) && (any(time_inst == k))
           x = z(1:n*(N+1));
           x = reshape(x,n,N+1);
           X1(iii,:) = x(1,:);
           X2(iii,:) = x(2,:);
           iii = iii + 1;
       end
        
       %% bar_z - update : Clipping on the box                           %%
  
       bar_z_old = bar_z;
       
       
       %%************************ Your code here ************************%%
       %%                                                                %%
       bar_z = min(max(z+mu,l),u);
       %%****************************************************************%%
        
       %% mu - update : Dual ascent step
       
       %%************************ Your code here ************************%%
       %%                                                                %%
       mu = mu + z - bar_z;
       %%****************************************************************%%
       
       %% diagnostics, reporting, termination checks    
       history.objval(k)  = 0.5 * z'*Mat*z;
       history.r_norm(k) = norm( z - bar_z);
       history.s_norm(k) = norm( rho * (bar_z - bar_z_old));
       history.eps_pri(k) = sqrt(n*(N+1)+m*N)*EPS_ABS + ...
                            EPS_REL*max( norm(bar_z), norm(z));
       history.eps_dual(k) = sqrt(n*(N+1)+m*N)*EPS_ABS + EPS_REL*norm(mu);
       if ~QUIET
            fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), history.s_norm(k),...
            history.eps_dual(k), history.objval(k));
       end
       
       if (history.r_norm(k)  < history.eps_pri(k) &&...
               history.s_norm(k) < history.eps_dual(k))%...
               % && k>max(time_inst))
             break;
       end   
    end
end