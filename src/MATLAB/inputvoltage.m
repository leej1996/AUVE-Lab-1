%% Solves Question 2.1
function u = inputvoltage(D,A,Delta,Ts)
% Creates linspace
t = 0:Ts:D;
% Creates square wave with amplitude A/2 and period D
u = A/2*square(2*pi*t/Delta);
end