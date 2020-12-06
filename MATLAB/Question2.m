%% Main function for the questions of the second part
function []= Question2()

%% Gets the input voltage
D = 500e-3; % duration
A = 0.1; % peak to peak amplitude
Delta = 100e-3; % period
Ts = 1e-3; % sampling time

u = inputvoltage(D,A,Delta,Ts);

%% Test
G = 50;
Ta = 20e-3;
Tf = 25e-3;
L = 512;

%% Symulations
%% When the system is perfectly modeled Ta=Tf
% Small q - > we trust our model
q = 0.0001;

% Initial state
x1_0 = [0.02;0];
P1_0 = [(2*pi/12)^2 0;0, 0];

[y,x] = simulate(u,G,Ta,Ts,L,x1_0);
xe = kal(y,u,G,Ta,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Ta,Ts,L,x1_0,q);

figure(1), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1)), plot(y);
title('theta');
figure(1), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)) %, plot(y);
title('omega');
suptitle('Perfect: KF with Small variance (q)');

figure(2), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1)), plot(y);
title('theta');
figure(2), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)) %, plot(y);
title('omega');
suptitle('Perfect: Stat KF with Small variance (q)');

% Big q -> we dont trust much our model
q = 0.2;

xe = kal(y,u,G,Ta,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Ta,Ts,L,x1_0,q);

figure(3), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1)), plot(y,'*');
title('theta');
figure(3), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)) %, plot(y);
title('omega');
suptitle('Perfect: KF with Big variance (q)');

figure(4), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1)), plot(y);
title('theta');
figure(4), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)) %, plot(y);
title('omega');
suptitle('Perfect: Stat KF with Big variance (q)');

%% When the model of the system is rough 
% Small q - > we trust our model
q = 0.0001;

[y,x] = simulate(u,G,Ta,Ts,L,x1_0);

xe = kal(y,u,G,Tf,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Tf,Ts,L,x1_0,q);

figure(5), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1)), plot(y);
title('theta');
figure(5), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)) %, plot(y);
title('omega');
suptitle('Rough: KF with Small variance (q)');

figure(6), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1)), plot(y);
title('theta');
figure(6), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)) %, plot(y);
title('omega');
suptitle('Rough: Stat KF with Small variance (q)');

% Big q -> we dont trust much our model
q = 0.2;

xe = kal(y,u,G,Tf,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Tf,Ts,L,x1_0,q);

figure(7), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1)), plot(y);
title('theta');
figure(7), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)) %, plot(y);
title('omega');
suptitle('Rough: KF with Big variance (q)');

figure(8), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1)), plot(y);
title('theta');
figure(8), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)) %, plot(y);
title('omega');
suptitle('Rough: Stat KF with Big variance (q)');

end

%% Solves Question 2.1
function u = inputvoltage(D,A,Delta,Ts)
% Creates linspace
t = 0:Ts:D;
% Creates square wave with amplitude A/2 and period D
u = A/2*square(2*pi*t/Delta);
end

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
