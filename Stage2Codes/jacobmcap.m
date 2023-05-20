function [J] = jacobmcap(q1,q2,lm1,lm2)
%JACOBMCAP Summary of this function goes here
%   Detailed explanation goes here
lm = [lm1 lm2];
J(1,1) = -lm(1)*sin(q1) - lm(2)*sin(q1+q2);
J(1,2) = -lm(2)*sin(q1+q2);
J(2,1) = lm(1)*cos(q1) + lm(2)*cos(q1+q2);
J(2,2) = -lm(2)*sin(q1+q2);
end

