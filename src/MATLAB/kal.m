%% Solves Question 2.3 - Kalman Filter
function xe = kal(y,u,G,T,Ts,L,xn,P,q)
% Matrices for the state-space continuos representation
A = [0,1;0,-1/T];
B = [0;G/T];
C = [1 0];
D = 0;

% Discretizing the system
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');

% For the Markov model
Hn = Cd;
Fn = Ad;

% Uniform distribution for the maesuremt
r = (2*pi/L/12)^2; %(a-b)/12

% Initializing the variable
xe = zeros(length(u),2);

% Loop of the algorithm 3.1
for i = 1:length(u)
    hn = Dd*u(i);
    fn = Bd*u(i);
    y_hat = Hn*xn + hn; % prediction
    Cxy = P*Hn';
    Cyy = Hn*P*Hn' + r;
    xn = xn + Cxy*inv(Cyy)*(y(i) - y_hat); % estimation
    xe(i,:) = xn';
    P = P - Cxy*inv(Cyy)*Cxy';
    xn = Fn*xn + fn;
    P = Fn*P*Fn' + Bd*q*Bd';
end
end