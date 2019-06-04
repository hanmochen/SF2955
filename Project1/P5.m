%% Load Data
load('HA1-data/stations.mat')
load('HA1-data/RSSI-measurements-unknown-sigma.mat')

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

logP=[];

States = [[0;0] [3.5;0] [0;3.5] [0;-3.5] [-3.5;0]];
Z=zeros(N,n+1);
for i=1:N
    Z(i,:) = simulate(mc,n);
end

%% SISR 
for v=0.5:0.05:3

    w = zeros(N,n);
    X = mvnrnd(zeros(6,1),diag([500,5,5,200,5,5]),N)';
    w0 = pdf(X,Y(:,1)',pos_vec,v);
    w(:,1) = w0;

    
    for  k = 1:(n-1) 
        % Update X
        X = pX* X + pZ*States(:,Z(:,k)) + pW*(mvnrnd([0,0],diag([0.25,0.25]),N)');
    
        %Resampling
        w0 =  pdf(X,Y(:,k+1)',pos_vec,v);  
        w(:,k+1) = w0;
    end

    p=sum(log(sum(w)));
    logP=[logP p];
end


%% Simulate with the maximum V   
    [~,ind]=max(logP);
    v=0.5:0.05:3;
    maxV=v(ind);
    figure,
    plot(0.5:0.05:3,logP),
    title('Log Probabilty - sigma')
    
    %% Initialization
    tau = zeros(2,n); 
    X = mvnrnd(zeros(6,1),diag([500,5,5,200,5,5]),N)';
    w0 = pdf(X,Y(:,1)',pos_vec,maxV);
    tau(1,1) = sum(X(1,:).*w0')/sum(w0);
    tau(2,1) = sum(X(4,:).*w0')/sum(w0);
    w(:,1) = w0;


    %% Main Loop
    for  k = 1:(n-1) 

       % Update X
        X = pX* X + pZ*States(:,Z(:,k)) + pW*(mvnrnd([0,0],diag([0.25,0.25]),N)');
    
        %Resampling
        w0 =  pdf(X,Y(:,k+1)',pos_vec,maxV);
        w(:,k+1) = w0;
        ind = randsample(N,N,true,w0);
        X = X(:,ind);
        
        %Update tau
        tau(1,k+1) = sum(X(1,:).*w0')/sum(w0);
        tau(2,k+1) = sum(X(4,:).*w0')/sum(w0);
    end

%% Draw Trajectory
figure,
plot(tau(1,:),tau(2,:));
hold on
plot(pos_vec(1,:),pos_vec(2,:),'*','Color','r');
title('Simulated Trajectory')

%% Calculate the observation PDF
function p=pdf(x,y,pos_vec,v)
    p=mvnpdf(y,90-30*log10(pdist2(x([1,4],:)',pos_vec')),v*v*diag(ones(6,1)));
end

