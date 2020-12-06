%% Solves Question 2.3 - Stationary Kalman Filter
function xe = stat_kal(y,u,G,T,Ts,L,x1_0,q)
% Matrices for the state-space continuos representation
A = [0,1;0,-1/T];
B = [0;G/T];
C = [1 0];
D = 0;

% Discretizing the system
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');

% Uniform distribution for the maesuremt
r = (2*pi/L/12)^2;
% K matrix
[K,~,~,~] = dlqe(Ad,Bd,Cd,q,r);

% For the Markov model
H = Cd;
F = Ad;

% Initializing the variables
xe = zeros(length(u),2);
xe(1,:) = x1_0';

% Loop of the algorithm 3.2
for i = 2:length(u)
    hn = Dd*u(i-1);
    fn = Bd*u(i-1);
    y_hat = H*xe(i-1,:)' + hn; % prediction
    x_hat = xe(i-1,:)' + K*(y(i-1) - y_hat); % estimation
    xe(i,:) = (F*x_hat + fn)';
end
end