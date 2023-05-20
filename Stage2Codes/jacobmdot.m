function Jdot = jacobmdot(q1,q2, q1dot, q2dot)
%JACOBMDOT Summary of this function goes here
%   Detailed explanation goes here
lm = [1.04 0.96];
Jdot(1,1) = -lm(1)*cos(q1)*q1dot-lm(2)*cos(q1+q2)*(q1dot + q2dot) ;
Jdot(1,2) = -lm(2)*cos(q1+q2)*(q1dot + q2dot);
Jdot(2,1) = -lm(1)*sin(q1)*q1dot-lm(2)*sin(q1+q2)*(q1dot + q2dot) ;
Jdot(2,2) = -lm(2)*sin(q1+q2)*(q1dot + q2dot);
end

