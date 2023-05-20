function Jdot = jacobmcapdot(q1,q2, q1dot, q2dot,lm1,lm2)
%JACOBMCAPDOT Summary of this function goes here
%   Detailed explanation goes here
lm = [lm1 lm2];
Jdot(1,1) = -lm(1)*cos(q1)*q1dot-lm(2)*cos(q1+q2)*(q1dot + q2dot) ;
Jdot(1,2) = -lm(2)*cos(q1+q2)*(q1dot + q2dot);
Jdot(2,1) = -lm(1)*sin(q1)*q1dot-lm(2)*sin(q1+q2)*(q1dot + q2dot) ;
Jdot(2,2) = -lm(2)*sin(q1+q2)*(q1dot + q2dot);
end

