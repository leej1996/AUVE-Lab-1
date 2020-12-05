function u = inputvoltage(D,A,Delta,Ts)
%Question 2.1
t = 0:Ts:D;
u = A/2*square(2*pi*t/Delta);
end

