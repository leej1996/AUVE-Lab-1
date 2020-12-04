function [y,x] = simulate(u,G,T,Ts,L,x1)
%Question 2.2
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

