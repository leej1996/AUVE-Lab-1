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
plot(x(:,1)), hold on,plot(xe(:,1)), plot(y); legend('Ground truth','Estimation','Measurments');
title('theta');
figure(1), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)); legend('Ground truth','Estimation');
title('omega');

sgtitle('Perfect: KF with Small variance (q)');

figure(2), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1)), plot(y); legend('Ground truth','Estimation','Measurments');
title('theta');
figure(2), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)); legend('Ground truth','Estimation');
title('omega');
sgtitle('Perfect: Stat KF with Small variance (q)');

% Big q -> we dont trust much our model
q = 0.005;

xe = kal(y,u,G,Ta,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Ta,Ts,L,x1_0,q);

figure(3), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1)), plot(y); legend('Ground truth','Estimation','Measurments');
title('theta');
figure(3), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)); legend('Ground truth','Estimation');
title('omega');
sgtitle('Perfect: KF with Big variance (q)');

figure(4), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1)), plot(y); legend('Ground truth','Estimation','Measurments');
title('theta');
figure(4), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)); legend('Ground truth','Estimation');
title('omega');
sgtitle('Perfect: Stat KF with Big variance (q)');

%% When the model of the system is rough 
% Small q - > we trust our model
q = 0.0001;

[y,x] = simulate(u,G,Ta,Ts,L,x1_0);

xe = kal(y,u,G,Tf,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Tf,Ts,L,x1_0,q);

figure(5), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1)), plot(y); legend('Ground truth','Estimation','Measurments');
title('theta');
figure(5), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)); legend('Ground truth','Estimation');
title('omega');
sgtitle('Rough: KF with Small variance (q)');

figure(6), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1)), plot(y); legend('Ground truth','Estimation','Measurments');
title('theta');
figure(6), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2)); legend('Ground truth','Estimation');
title('omega');
sgtitle('Rough: Stat KF with Small variance (q)');

% Big q -> we dont trust much our model
q = 0.005;

xe = kal(y,u,G,Tf,Ts,L,x1_0,P1_0,q);
xe_s = stat_kal(y,u,G,Tf,Ts,L,x1_0,q);

figure(7), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe(:,1)), plot(y); legend('Ground truth','Estimation','Measurments');
title('theta');
figure(7), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe(:,2)); legend('Ground truth','Estimation');
title('omega');
sgtitle('Rough: KF with Big variance (q)');

figure(8), subplot(211), hold off
plot(x(:,1)), hold on,plot(xe_s(:,1)), plot(y); legend('Ground truth','Estimation','Measurments');
title('theta');
figure(8), subplot(212), hold off
plot(x(:,2)), hold on, plot(xe_s(:,2));  legend('Ground truth','Estimation');
title('omega');
sgtitle('Rough: Stat KF with Big variance (q)');

end