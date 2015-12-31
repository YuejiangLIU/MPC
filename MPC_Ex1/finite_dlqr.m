function [K0,V0] = finite_dlqr(A,B,Q,R,Pf,N)
% MPC Ex1 | Yuejiang LIU | Mar. 2015

K = zeros(size(B));
H = Pf;

for i = N-1:-1:0
    K_new = - inv(R+B'*H*B)*B'*H*A;
    H_new = Q + K_new'*R*K_new + (A+B*K_new)'*H*(A+B*K_new);
    K = K_new;
    H = H_new;
end

K0 = K;
V0 = H; 

end
