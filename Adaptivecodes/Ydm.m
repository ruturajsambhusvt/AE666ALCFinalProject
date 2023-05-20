function [Y] = Ydm(q1, q2,q1dot, q2dot, qr1dot, qr2dot, qr1ddot, qr2ddot)
%YDM Summary of this function goes here
%   Detailed explanation goes here
Y(1,1) = qr1ddot;
Y(1,2) = 2*cos(q2)*qr1ddot + cos(q2)*qr2ddot - sin(q2)*q2dot*qr1dot - sin(q2)*(q1dot + q2dot)*qr2dot;
Y(1,3) = qr2ddot;
Y(2,1)  = 0;
Y(2,2) = cos(q2)*qr1ddot + sin(q2)*q1dot*qr1dot;
Y(2,3) = qr1ddot + qr2ddot;
end

