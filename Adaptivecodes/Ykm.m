function [Y] = Ykm(q1, q2,q1dot, q2dot)
%YKM Summary of this function goes here
%   Detailed explanation goes here
Y(1,1) = -sin(q1)*q1dot;
Y(1,2) = -sin(q1+q2)*(q1dot+q2dot);
Y(2,1) = cos(q1)*q1dot;
Y(2,2) = cos(q1+q2)*(q1dot+q2dot);
end

