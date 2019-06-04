%% Load Data
clc;clear all;close all;
load('HA1-data/RSSI-measurements.mat')
load('HA1-data/stations.mat')

%% Basic Settings

%Constants
N = 10000;
n = 501;
dT = 0.5;
alpha = 0.6;

% Transition probability matrix
P = 1/20*(15*diag(ones(1,5))+ones(5));

% Matrices 
phiX = [1 dT (dT^2)/2;...
        0 1 dT;...
        0 0 alpha];
pX = [phiX zeros(3,3);...
       zeros(3,3) phiX];
phiZ = [(dT^2)/2;...
        dT;...
        0];
pZ = [phiZ zeros(3,1); ...
    zeros(3,1) phiZ];
phiW = [(dT^2)/2;...
        dT;...
        1];
pW = [phiW zeros(3,1);...
    zeros(3,1) phiW];


%% Initialization
tau = zeros(2,n); 
w = zeros(N,n);

X = mvnrnd(zeros(6,1),diag([500,5,5,200,5,5]),N)';
w0 = pdf(X,Y(:,1)',pos_vec);
tau(1,1) = sum(X(1,:).*w0')/sum(w0);
tau(2,1) = sum(X(4,:).*w0')/sum(w0);
w(:,1) = w0;

States = [[0;0] [3.5;0] [0;3.5] [0;-3.5] [-3.5;0]];
mc = dtmc(P);
Z=zeros(N,n);
for i=1:N
    Z(i,:) = simulate(mc,n-1);
end

%% Loop
for  k = 1:(n-1) 
    %Update X
    X = pX* X + pZ*States(:,Z(:,k)) + pW*(mvnrnd([0,0],diag([0.25,0.25]),N)');
    
    %Update Weights
    w0 =  w0.*pdf(X,Y(:,k+1)',pos_vec);
    w0 = w0/sum(w0);
    w(:,k+1) = w0;
    
    %Estimate tau
    tau(1,k+1) = sum(X(1,:).*w0')/sum(w0);
    tau(2,k+1) = sum(X(4,:).*w0')/sum(w0);

end

%% Draw Trajectory
figure,
plot(tau(1,:),tau(2,:),'.');
hold on
plot(pos_vec(1,:),pos_vec(2,:),'*','Color','r');
title('Simulated Trajectory')


%% Calculating the efficency sampling
track_sample = zeros(1,n);
sample_size = 0;
dummy = Inf;
for i = 1:n
    CV2 = (1/N)*sum((N*(w(:,i)./sum(w(:,i)))-1).^2);
    track_sample(i) = N/(1+CV2);
    if track_sample(i)<dummy
        dummy = track_sample(i);
        sample_size = i;       
    end
end

%% Plot histogram

%n= 50
figure,subplot(3,1,1),
histogram(log(w(:,50)),20);
title('n=50')
%n= 100
subplot(3,1,2),
histogram(log(w(:,100)),20);
title('n=100')
%n= 200
subplot(3,1,3),
histogram(log(w(:,200)),20);
title('n=200')

%% Plot efficient sample size
figure,
plot(track_sample)

%% Calculate the observation PDF
function p=pdf(x,y,pos_vec)
    p=mvnpdf(y,90-30*log10(pdist2(x([1,4],:)',pos_vec')),diag(2.25*ones(6,1)));
end

