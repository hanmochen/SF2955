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

prob=[];
%% SISR
for v=0.5:0.05:3

%% Initialization
    tau = zeros(2,n); 
    tracK = zeros(N,n);

    X = mvnrnd(zeros(6,1),diag([500,5,5,200,5,5]),N)';
    w = pdf(X,Y(:,1)',pos_vec,v);
    tau(1,1) = sum(X(1,:).*w')/sum(w);
    tau(2,1) = sum(X(4,:).*w')/sum(w);
    track(:,1) = w;

    States = [[0;0] [3.5;0] [0;3.5] [0;-3.5] [-3.5;0]];
    mc = dtmc(P);
    simulate_Z = simulate(mc,n);

    %% Main Loop
    for  k = 1:(n-1) 
        zM = repmat(pZ*States(:,simulate_Z(k)),1,N);
        wM = pW*(mvnrnd([0,0],diag([0.25,0.25]),N)');
        xM = pX* X; 
        YmN = repmat(Y(:,k+1),1,N);
        X = xM + zM + wM;
        w =  pdf(X,Y(:,k+1)',pos_vec,v);
        if w==0
           p=-Inf;
           prob=[prob p];
           break
        else
            ind = randsample(N,N,true,w);
            X = X(:,ind);
        end 
        tau(1,k+1) = sum(X(1,:).*w')/sum(w);
        tau(2,k+1) = sum(X(4,:).*w')/sum(w);
        track(:,k+1) = w;
    end
    
    Z=mvnpdf(Y',90-30*log10(pdist2(tau',pos_vec')),v*v*diag(ones(6,1)));
    p=sum(log(Z));
    prob=[prob p];
end

%% Simulate with the maximum V   
   [~,ind]=max(prob);
    v=0.5:0.05:3;
    maxV=v(ind);

    %% Initialization
    tau = zeros(2,n); 
    tracK = zeros(N,n);

    X = mvnrnd(zeros(6,1),diag([500,5,5,200,5,5]),N)';
    w = pdf(X,Y(:,1)',pos_vec,maxV);
    tau(1,1) = sum(X(1,:).*w')/sum(w);
    tau(2,1) = sum(X(4,:).*w')/sum(w);
    track(:,1) = w;

    States = [[0;0] [3.5;0] [0;3.5] [0;-3.5] [-3.5;0]];
    mc = dtmc(P);
    simulate_Z = simulate(mc,n);

    %% Main Loop
    for  k = 1:(n-1) 
        zM = repmat(pZ*States(:,simulate_Z(k)),1,N);
        wM = pW*(mvnrnd([0,0],diag([0.25,0.25]),N)');
        xM = pX* X; 
        YmN = repmat(Y(:,k+1),1,N);
        X = xM + zM + wM;
        w =  pdf(X,Y(:,k+1)',pos_vec,maxV);
        ind = randsample(N,N,true,w);
        X = X(:,ind);
        tau(1,k+1) = sum(X(1,:).*w')/sum(w);
        tau(2,k+1) = sum(X(4,:).*w')/sum(w);
        track(:,k+1) = w;
    end

%% Draw Trajectory
figure,
plot(tau(1,:),tau(2,:),'LineWidth',2);
hold on
plot(pos_vec(1,:),pos_vec(2,:),'*','Color','r');
title('Simulated Trajectory')

%% Calculate the observation PDF
function p=pdf(x,y,pos_vec,v)
    p=mvnpdf(y,90-30*log10(pdist2(x([1,4],:)',pos_vec')),v*v*diag(ones(6,1)));
end

