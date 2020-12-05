function []= Question2()




G = 50;
Ta = 20e-3;
Tf = 25e-3;
Ts = 1e-3;
L = 512;

D = 500e-3;
A = 0.1;
Delta = 100e-3;

%% System is Perfect
%small q
q = 0.000000000001;
x1_0 = [0.02;0];
P1_0 = [(2*pi/12)^2 0;0, 0];

u = inputvoltage(D,A,Delta,Ts);
[y,x] = simulate(u,G,Ta,Ts,L,x1_0);

xe = kal(y,u,G,Ta,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Ta,Ts,L,x1_0,q);

figure(1), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1))%, plot(y);
title('theta');
figure(1), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)) %, plot(y);
title('omega');
sgtitle('Perfect: KF with Small variance (q)');

figure(2), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1))%, plot(y);
title('theta');
figure(2), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)) %, plot(y);
title('omega');
sgtitle('Perfect: Stat KF with Small variance (q)');

%Big q
q = 0.00000001;

xe = kal(y,u,G,Ta,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Ta,Ts,L,x1_0,q);

figure(3), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1))%, plot(y);
title('theta');
figure(3), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)) %, plot(y);
title('omega');
sgtitle('Perfect: KF with Big variance (q)');

figure(4), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1))%, plot(y);
title('theta');
figure(4), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)) %, plot(y);
title('omega');
sgtitle('Perfect: Stat KF with Big variance (q)');

%% System is Rough
%small q
q = 0.000000000001;
x1_0 = [0.02;0];
P1_0 = [(2*pi/12)^2 0;0, 0];

u = inputvoltage(D,A,Delta,Ts);
[y,x] = simulate(u,G,Ta,Ts,L,x1_0);

xe = kal(y,u,G,Tf,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Tf,Ts,L,x1_0,q);

figure(5), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1))%, plot(y);
title('theta');
figure(5), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)) %, plot(y);
title('omega');
sgtitle('Rough: KF with Small variance (q)');

figure(6), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1))%, plot(y);
title('theta');
figure(6), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)) %, plot(y);
title('omega');
sgtitle('Rough: Stat KF with Small variance (q)');

%Big q
q = 0.00000001;

xe = kal(y,u,G,Tf,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Tf,Ts,L,x1_0,q);

figure(7), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1))%, plot(y);
title('theta');
figure(7), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)) %, plot(y);
title('omega');
sgtitle('Rough: KF with Big variance (q)');

figure(8), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1))%, plot(y);
title('theta');
figure(8), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)) %, plot(y);
title('omega');
sgtitle('Rough: Stat KF with Big variance (q)');

end

%% Solves Question 2.1
function u = inputvoltage(D,A,Delta,Ts)
t = 0:Ts:D;
u = A/2*square(2*pi*t/Delta);
end

%% Solves Question 2.2
function [y,x] = simulate(u,G,T,Ts,L,x1)
A = [0,1;0,-1/T];
B = [0;G/T];
C = [1 0];
D = 0;
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');
x = zeros(length(u),2);
theta = zeros(1,length(u));
x(1,:) = x1';
theta(1) = x1(1);
for i = 2:length(u)
    x(i,:) = (Ad*x(i-1,:)' + Bd*u(i-1))';
    theta(i) = Cd*x(i-1,:)' + Dd*u(i-1);
end
y = round(theta.*L/2/pi)*2*pi/L;
end

function xe = kal(y,u,G,T,Ts,L,xn,P,q)
% Question 2.3
A = [0,1;0,-1/T];
B = [0;G/T];
C = [1 0];
D = 0;
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');
Hn = Cd;
Fn = Ad;

r = (2*pi/L/12)^2; %(a-b)/12
xe = zeros(length(u),2);

for i = 1:length(u)
    hn = Dd*u(i);
    fn = Bd*u(i);
    y_hat = Hn*xn + hn;
    Cxy = P*Hn';
    Cyy = Hn*P*Hn' + r;
    xn = xn + Cxy*inv(Cyy)*(y(i) - y_hat);
    xe(i,:) = xn';
    P = P - Cxy*inv(Cyy)*Cxy';
    xn = Fn*xn + fn;
    P = Fn*P*Fn' + B*q*B';
end
end

function xe = stat_kal(y,u,G,T,Ts,L,x1_0,q)
% Question 2.3
A = [0,1;0,-1/T];
B = [0;G/T];
C = [1 0];
D = 0;
r = (2*pi/L/12)^2;
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');
[K,~,~,~] = dlqe(Ad,Bd,Cd,q,r);

H = Cd;
F = Ad;

xe = zeros(length(u),2);
xe(1,:) = x1_0';

for i = 2:length(u)
    hn = Dd*u(i-1);
    fn = Bd*u(i-1);
    y_hat = H*xe(i-1,:)' + hn;
    x_hat = xe(i-1,:)' + K*(y(i-1) - y_hat);
    xe(i,:) = (F*x_hat + fn)';
end
end


