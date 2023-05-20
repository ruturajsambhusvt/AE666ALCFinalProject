function [Mm] = inertias(q1,q2)
%INERTIAM Summary of this function goes here
%   Detailed explanation goes here
%robot2
massm = [2.54 1.82];
Im = [0.12 0.06];
lm = [1.28 0.84];
Mm(1,1) = massm(1)*(lm(1)/2)^2 + massm(2)*(lm(1)^2 + (lm(2)/2)^2 + 2*lm(1)*(lm(2)/2)*cos(q2)) + Im(1)  + Im(2);
Mm(1,2) = massm(2)*((lm(2)/2)^2 + lm(1)*(lm(2)/2)*cos(q2)) + Im(2);
Mm(2,1) = massm(2)*((lm(2)/2)^2 + lm(1)*(lm(2)/2)*cos(q2)) + Im(2);
Mm(2,2) = massm(2)*(lm(2)/2)^2 + Im(2);
end