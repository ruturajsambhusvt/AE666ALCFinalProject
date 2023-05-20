function [J] = jacobm(q1,q2)
%JACOB Summary of this function goes here
%   Detailed explanation goes here
lm = [1.04 0.96];
J(1,1) = -lm(1)*sin(q1) - lm(2)*sin(q1+q2);
J(1,2) = -lm(2)*sin(q1+q2);
J(2,1) = lm(1)*cos(q1) + lm(2)*cos(q1+q2);
J(2,2) = -lm(2)*sin(q1+q2);
end

