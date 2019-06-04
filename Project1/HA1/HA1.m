%%Home Assignment 1
%% Problem 1 - Moving Model
clear all; close all; clc;
%Trajecory lenght
m = 200;

%Initil values for Z
Zvalues = [0 , 0; 3.5, 0; 0, 3.5; 0, -3.5; -3.5, 0];
Zindex = randi(5,1);
Z0 = Zvalues(Zindex,:);

%Transition Probabilities of Z
P = (1/20)*[16 1 1 1 1; 1 16 1 1 1; 1 1 16 1 1; 1 1 1 16 1; 1 1 1 1 16];
 
%Initial values for X
mu=zeros(6,1);
sig=diag([500,5,5,200,5,5]);
X0 = mvnrnd(mu, sig);

%Values for W
sigma1= 0.5;
sigW = diag([sigma1^2 sigma1^2]);
muW = [0;0];
W = mvnrnd(muW, sigW, m)';

%Values of PHI, PSI_z and PSI_w
dt = 0.5; alpha=0.6;
phiT = [1 dt (dt^2)/2; 0 1 dt; 0 0 alpha];
psiTZ = [(dt^2)/2; dt; 0]; psiTW = [(dt^2)/2; dt; 1];

PHI = [phiT zeros(3,3); zeros(3,3) phiT];
PSI_z = [psiTZ zeros(3,1); zeros(3,1) psiTZ];
PSI_w = [psiTW zeros(3,1); zeros(3,1) psiTW];

%Trajectory calulation
n = 1000;
X = zeros(6,m);
X(:,1) = X0;
Z = Z0';
Zprob = rand(1,m);

for j=2:m
    X(:,j) = PHI*X(:,j-1) + PSI_z*Z + PSI_w*W(:,j);
    Zindex = Zfunc(Zindex, P, Zprob);
    Z = Zvalues(Zindex,:)';  
end 

Trajectory_cor = [X(1,:); X(4,:)];
plot(Trajectory_cor(1,:), Trajectory_cor(2,:))
 
%% probelm 3

clc;
load('stations.mat') %pos_vec is 6 positions for the basis stations.
load('RSSI-measurements.mat') % Y is the RSSI that at mobile unit receive from the lth BS at time n 

%Creating 10000 random samples for the initial state vector X0
mu=zeros(6,1); sig=diag([500,5,5,200,5,5]); N = 10000; 
X = mvnrnd(mu,sig,N)';
Xcor = [X(1,:);X(4,:)]; %coordinates for trajectory

%parameters
varsig = 1.5;
sigma2 = eye(6)*varsig^2; 
m = 500;

% Weights
initialweights = mvnpdf(Y(:,1)', MUfunc(pos_vec, Xcor, N), sigma2);
w = zeros(m, N);
w(1,:) = initialweights;

% tau_1=t1 & tau_2=t2, initial value for tau = t0;
t1 = zeros(1,m);
t2 = zeros(1,m);
t0 = sum(Xcor.*repmat(initialweights',2,1),2)/sum(initialweights'); %initial value for tau
Tau = [t1; t2];
Tau(:,1) = t0;
 
% Values for Z
Zvalues = [0 , 0; 3.5, 0; 0, 3.5; 0, -3.5; -3.5, 0];
Zindex = zeros(N,m);
Zindex(:,1) = randi(5,1,N);
Zvalues = Zvalues(Zindex(:,1),:)';
Z = Zvalues(:,1);

%SIS algorithm
for i = 2:m
    %new noise
    W = mvnrnd([0;0], eye(2)*sigma1^2,N)';
    
    %new X coordinates
    X = PHI*X + PSI_z*Z + PSI_w*W; 
    Xcor = [X(1,:);X(4,:)];

    %new weight
    w(i,:) = w(i-1,:).*(mvnpdf(Y(:,i)', MUfunc(pos_vec, Xcor, N),sigma2))'; 
    
    %calculating estimated position Tau
    Tau(:,i) = sum(Xcor.*repmat(w(i,:),2,1),2)/sum(w(i,:)); 
    
    %new Z
    Zprob = rand(1,N); 
    Zindex(:,i) = Zfunc2(Zindex(:,i-1),P,Zprob); 
    Z = Zvalues(:,Zindex(:,i));
end

% Plots Problem 3

%Plot of only the path
figure
plot(Tau(1,:),Tau(2,:))
xlabel('X1')
ylabel('X2')

%Plot of path and basis staions
figure
plot(Tau(1,:),Tau(2,:))
hold on
plot(pos_vec(1,:),pos_vec(2,:),'ob')
xlabel('X1')
ylabel('X2')

%Histogram of weights
weights = [1 5 10];
figure
for i = 1:length(weights)
    subplot(3,1,i)
    histogram(log10(w(weights(i),:)));    
    xlabel(['weight{', num2str(weights(i)),'}']) 
end

% Efficient sample size
Tstep = 1:m;
ESS = zeros(1, length(Tstep));
for i=1:length(Tstep)
    CV = (1/N) * sum(((N.*w(Tstep(i),:)./sum(w(Tstep(i),:))) - ones(1, N)).^2);
    ESS(i) = N / (1 + (CV.^2));
end

% Plot of efficient sample size
figure
plot(Tstep, ESS)
xlabel('Time')
ylabel('Efficient sample size')

%% Problem 4
clc;
%initial values of X
N=10000;
X = mvnrnd(mu,sig,N)';
Xcor = [X(1,:);X(4,:)];

% Weights
initialweights = mvnpdf(Y(:,1)', MUfunc(pos_vec, Xcor, N), sigma2);
w = zeros(m, N);
w(1,:) = initialweights;

% tau_1=t1 & tau_2=t2, initial value for tau = t0;
t1 = zeros(1,m);
t2 = zeros(1,m);
t0 = sum(Xcor.*repmat(initialweights',2,1),2)/sum(initialweights'); %initial value for tau
Tau = [t1; t2];
Tau(:,1) = t0; 

% Zindex = zeros(N,m);
% Zindex(:,1) = randi(5,1,N); 
Zindex = randi(5,1);
Zvalues = [0 , 0; 3.5, 0;0, 3.5;0, -3.5;-3.5, 0];      
Z = Zvalues(Zindex,:);
Z = repmat(Z,2,N/2);

MostProbZ = ones(m-1,2);

for i = 2:m
    %new noise
    W = mvnrnd([0;0], eye(2)*sigma1^2,N)';
  
    %new X coordinates
    X = PHI*X + PSI_z*Z + PSI_w*W;
    Xcor = [X(1,:);X(4,:)];
    
    %calculate weights
    w(i,:) =  (mvnpdf(Y(:,i)', MUfunc(pos_vec,Xcor,N), sigma2))'; 
    
    % Estimation of Tau
    Tau(:,i) = sum(Xcor.*repmat(w(i,:),2,1),2)/sum(w(i,:)); 
    ind = randsample(N,N,true,w(i,:)); % selection
    X = X(:,ind);
    
    %New Z
    Zprob = rand(1);
    Zindex = Zfunc3(Zindex,P,Zprob);
    Z = Zvalues(Zindex,:);
    Z = repmat(Z,2,N/2);
       
    %most probable driving command
    Zprob2 = P(1,1)-0.01;
    Zindex2 = Zindex;
    Zindex2 = Zfunc3(Zindex2,P,Zprob2);
    MostProbZ(i-1,:) = Zvalues(Zindex2,:)'; 
end

% Plots Problem 4

%only the path
figure
plot(Tau(1,:),Tau(2,:))
xlabel('X1')
ylabel('X2')

%path with basis staions
figure
plot(Tau(1,:),Tau(2,:))
hold on
plot(pos_vec(1,:),pos_vec(2,:),'ob')
xlabel('X1')
ylabel('X2')


% %Histogram of weights
% weights = [1 5 10];
% figure
% for i = 1:length(weights)
%     subplot(3,1,i)
%     histogram(log10(w(weights(i),:)));    
%     xlabel(['weight{', num2str(weights(i)),'}']) 
% end
 
% % Efficient sample size
% Tstep = 1:m;
% ESS = zeros(1, length(Tstep));
% for i=1:length(Tstep)
%     CV = (1/N) * sum(((N.*w(Tstep(i),:)./sum(w(Tstep(i),:))) - ones(1, N)).^2);
%     ESS(i) = N / (1 + (CV.^2));
% end
% 
% % Plot of efficient sample size
% figure
% plot(Tstep, ESS)
% xlabel('Time')
% ylabel('Efficient sample size')

%% Problem 5
close all; clc;
load('RSSI-measurements-unknown-sigma.mat')
load('stations.mat') 

vs = 0.1:0.1:3; %generated grid for varsigma
k = length(vs); Tau = cell(1,k); N = 10000;
logLH = zeros(k,1); %log-likelihood

for j= 1:k
    
    %initial value for X
    X = mvnrnd(mu,sig,N)';
    Xcor = [X(1,:);X(4,:)];
    
    %initial values for Z
    Zvalues = [0 , 0; 3.5, 0; 0, 3.5; 0, -3.5; -3.5, 0];
    Zindex = zeros(N,m);
    Zindex(:,1) = randi(5,1,N);
    Zvalues = Zvalues(Zindex(:,1),:)';
    Z = Zvalues(:,1);
    
    % Weights
    initialweights = mvnpdf(Y(:,1)', MUfunc(pos_vec, Xcor, N), eye(6)*vs(j));
    w = zeros(m, N);
    w(1,:) = initialweights;

    % t1 & t2, initial value for tau = t0;
    t1 = zeros(1,m); t2 = zeros(1,m);
    t0 = sum(Xcor.*repmat(initialweights',2,1),2)/sum(initialweights');
    Tau{j} = [t1; t2];
    Tau{j}(:,1) = t0; 
   
    for i = 2:m
        %new noise
        W = mvnrnd([0;0], eye(2)*sigma1^2,N)';

        %new X coordinates
        X = PHI*X + PSI_z*Z + PSI_w*W;
        Xcor = [X(1,:);X(4,:)];

        %calculate weights
        w(i,:) =  (mvnpdf(Y(:,i)', MUfunc(pos_vec,Xcor,N), eye(6)*vs(j)))'; 

        % Estimation of Tau
        Tau{j}(:,i) = sum(Xcor.*repmat(w(i,:),2,1),2)/sum(w(i,:)); 
        
        %This if w is non-negative or zero
        w2 = abs(w(i,:));
        if w2 ~= zeros(1,length(w2(1,:)))
            ind = randsample(N,N,true,w2); % selection
            X = X(:,ind);
        end

        %new Z
        Zprob = rand(1,N); 
        Zindex(:,i) = Zfunc2(Zindex(:,i-1),P,Zprob); 
        Z = Zvalues(:,Zindex(:,i));
    end
    
    C = 0;
    for i = 1:m
    C = C + mean(log((1/N)*sum(w(i,:))));    
    end
    logLH(j) = 1/m*C;
end


% PLot of Tau for the estimated varsigma
vs_hat = vs(find(logLH==max(logLH)));
disp(['Best estimation of varsigma is ' num2str(vs_hat)]);
i = find(logLH==max(logLH));
plot(Tau{i}(1,:),Tau{i}(2,:),'-')
hold on
plot(pos_vec(1,:),pos_vec(2,:),'ob')


