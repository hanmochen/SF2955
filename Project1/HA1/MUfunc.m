function [ mu ] = MUfunc(pos_vec, X, N)

v = 90; eta = 3;

%Replicate to avoid for-loops
X0cor = repmat(X,6,1); % replicate the X0 coordinates 
BScor = repmat(pos_vec(:),1,N); %replication of the stations positions

d = (X0cor - BScor)';
eucDist = zeros(N,6);

for i = 1:6
       
       eucDist(:,i) = sqrt(abs(d(:,2 * i-1)).^2 + abs(d(:,2 * i)).^2); %Euclidean distance
end

mu = v - 10*eta*log10(eucDist);

end
