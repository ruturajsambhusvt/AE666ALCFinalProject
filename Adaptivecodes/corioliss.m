function [Cm] = corioliss(q1,q2,q1dot,q2dot)
%CORIOLISM Summary of this function goes here
%   Detailed explanation goes here
%robot2
massm = [2.54 1.82];
Im = [0.12 0.06];
lm = [1.28 0.84];
h = -massm(2)*lm(1)*(lm(2)/2)*sin(q2);
Cm(1,1) = h*q2dot;
Cm(1,2) = h*q2dot + h*q1dot;
Cm(2,1) = -h*q1dot;
Cm(2,2) = 0;
end


