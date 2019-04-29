m = 1000;%m steps
%% Basic Settings 

dT = 0.5;
alpha = 0.6;
X0 = (mvnrnd(zeros(6,1),diag([500,5,5,200,5,5])))';

% States
Z_state = [[0;0] [3.5;0] [0;3.5] [0;-3.5] [-3.5;0]];

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


%% Simulate a trajectory
% Simulate Z
mc = dtmc(P);
Z = simulate(mc,m);%simulate Z from a random state for m steps

% Simulate X
X = zeros(6,m+1);
X(:,1) = X0;
for i=1:m
    X(:,i+1) = pX* X(:,i) + pZ*Z_state(:,Z(i)) + pW*(mvnrnd(zeros(2,1),diag([0.25,0.25])))';
end

%% Draw the trajectory

figure
l= animatedline('Color','r','LineWidth',1);
minX = min(X(1,:));
minY = min(X(4,:));
maxX = max(X(1,:));
maxY = max(X(4,:));
axis([minX maxX minY maxY])
for i=1:(m+1)
    addpoints(l,X(1,i),X(4,i))
    drawnow    
end
