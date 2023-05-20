function [Mm] = inertiam(q1,q2)
%INERTIAM Summary of this function goes here
%   Detailed explanation goes here
%robot1
massm = [3.14 2.26];
Im = [0.16 0.07];
lm = [1.04 0.96];
Mm(1,1) = massm(1)*(lm(1)/2)^2 + massm(2)*(lm(1)^2 + (lm(2)/2)^2 + 2*lm(1)*(lm(2)/2)*cos(q2)) + Im(1)  + Im(2);
Mm(1,2) = massm(2)*((lm(2)/2)^2 + lm(1)*(lm(2)/2)*cos(q2)) + Im(2);
Mm(2,1) = massm(2)*((lm(2)/2)^2 + lm(1)*(lm(2)/2)*cos(q2)) + Im(2);
Mm(2,2) = massm(2)*(lm(2)/2)^2 + Im(2);
end

