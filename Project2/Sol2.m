%% Load Data
clear;clc;close all;
y = csvread("HA2-data/mixture-observations.csv");

%% Initialize
n = 1000;
N = 100;
theta = zeros(N,1);
theta0 = 0.5;% or  theta0 = rand

%% EM algorithm
for j = 1:N
    % E step
    p1 = normpdf(y,1,2);
    p0 = normpdf(y,0,1);
    p = (theta0*p1)./(theta0*p1 +(1-theta0)*p0);
    
    % M step
    theta0 = 1/n*(sum(p));
    theta(j) = theta0;
end

%% Plot Estimator
figure
plot(1:N,theta); 
title('$\hat\theta$','Interpreter','latex')
xlabel('Iteration times')
ylabel('$\hat\theta$','Interpreter','latex')

disp("The estimated theta is "+num2str(mean(theta(80:100))))