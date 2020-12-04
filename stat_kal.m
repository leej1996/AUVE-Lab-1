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

