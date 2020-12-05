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

