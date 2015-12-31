% 
clear all;
close all;
clc;

N_point = 20;
mu_series = logspace(-7,-1,N_point); 

N_repeat = 5;
    
avg_mu = zeros(N_point,1);
avg_gap = zeros(N_point,1);
avg_iter = zeros(N_point,1);
std_gap = zeros(N_point,1);
std_iter = zeros(N_point,1);

for i = 1:N_point 
    mu = zeros(N_repeat,1);
    gap = zeros(N_repeat,1);
    iter = zeros(N_repeat,1);
    for repeat = 1:N_repeat;
        result = ex2(mu_series(i));
        mu(repeat) = result{1,1};
        gap(repeat) = result{1,2};
        iter(repeat) = result{1,3};
    end
    avg_mu(i) = mean(mu);
    avg_gap(i) = mean(gap);
    avg_iter(i) = mean(iter);
    std_gap(i) = std(gap);
    std_iter(i) = std(iter);
end

avg = table(mu_series',avg_mu,avg_gap,avg_iter,std_gap,std_iter)

figure; 
errorbar(avg_mu,avg_iter,std_iter,'r');
set(gca,'xscale','log');
xlabel('\mu');
ylabel('Newton Iterations');
title('EPFL MPC Ex2 | Yuejiang LIU | March 2015');
grid on;
legend('Dim = 2');

writetable(avg,'dim2.dat');