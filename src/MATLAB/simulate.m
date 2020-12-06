%% Solves Question 2.2
function [y,x] = simulate(u,G,T,Ts,L,x1)
% Matrices for the state-space continuos representation
A = [0,1;0,-1/T];
B = [0;G/T];
C = [1 0];
D = 0;

% Discretizing the system
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');

%% Simulating the discrete system
% Initializing variables
x = zeros(length(u),2);  
theta = zeros(length(u),1);

% Assigning initial values 
x(1,:) = x1'; % also theta
theta(1) = x1(1);

% Computing the state variables and output
for i = 2:length(u)
    x(i,:) = Ad*x(i-1,:)' + Bd*u(i-1);
    theta(i) = Cd*x(i-1,:)' + Dd*u(i-1);
end

y = round(theta.*L/2/pi)*2*pi/L;
end