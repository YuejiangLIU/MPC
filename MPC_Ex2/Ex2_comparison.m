% load 
close all;
clear all;
clc;

data2 = importdata('dim2.dat');
data50 = importdata('dim50.dat');

figure; 
plot1 = errorbar(data2.data(:,1),data2.data(:,4),data2.data(:,6),'b');
hold on;
plot2 = errorbar(data50.data(:,1),data50.data(:,4),data50.data(:,6),'r');
set(gca,'xscale','log');
xlabel('\mu');
ylabel('Newton Iterations');
title('EPFL MPC Ex2 | Yuejiang LIU | March 2015');
grid on;
legend(([plot1, plot2]),'Dim = 2','Dim = 50');

