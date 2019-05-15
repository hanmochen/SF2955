%% Load Data
clear;clc;close all;
T = csvread("HA2-data/coal-mine.csv");

%% Initialize
%hyperparameter
var_theta = 1;
rho = 1;

d = 3;%number of break points
M = 10000;%steps

t = zeros(d+1,M);%t1,t2,...,t(d+1)
t(1,:) = 1851;
t(d+1,:) = 1963;
t(:,1)=1851:(1963-1851)/d:1963; % Evenly divided by d

theta = zeros(1,M);%theta
lambda = zeros(d,M);%lambda
theta(1) = gamrnd(2,var_theta);%prior
lambda(:,1) = gamrnd(2,theta(1),d,1);%prior


%% hybrid MCMC
for j=1:M-1
    
    %update theta
    theta(j+1) = gamrnd(2d+2,1/(sum(lambda(:,j)) + var_theta));
    
    %update lambda
    nT=histcounts(T,t(:,j));%Calculate n(t)
    lambda(:,j+1)=gamrnd(2+nT(:),1./(diff(t(:,j))+ theta(j+1)));
    
    %update t using Random walk proposal
    for i=2:d
        R = rho*(t(i+1,j)-t(i-1,j+1));
        temp = t(i,j)+ 2*R*rand -R;
        if temp > t(i+1,j) || temp < t(i-1,j+1) % check the low and high bound
            alpha = 0;
        else
            alpha = min(1,((temp-t(i-1,j+1))*(t(i+1,j)-temp))/...
                ((t(i,j) - (t(i-1,j+1)))*(t(i+1,j) - (t(i,j)))));%(t2'-t1)(t3-t2')/(t2-t1)(t3-t2)
        end
        
        if rand<=alpha
            t(i,j+1) = temp;
        else
            t(i,j+1)=t(i,j);
        end
    end    
end
%% Figure

figure(1);
for i=2:d
    histogram(t(i,:));
    title("Histogram for t")
    legend('t2','t3')
    hold on;
    xlabel("t")
    ylabel("Frequency")
end


figure(2);
for i=1:d
    histogram(lambda(i,:),0:0.1:2);
    hold on;
    title("Histogram for $\lambda$ ",'Interpreter','latex')
    xlabel("$\lambda$",'Interpreter','latex')
    ylabel("Frequency")        
end

figure(3);
histogram(theta,20);
title("Histogram for $\theta$ ",'Interpreter','latex')
xlabel("$\theta$",'Interpreter','latex')
ylabel("Frequency")   
    
