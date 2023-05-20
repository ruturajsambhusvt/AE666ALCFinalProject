function [Dyn] = UnconRobot(t,X)
%UNCONTROLLEDROBOT Summary of this function goes here
%   Detailed explanation goes here

q1m = X(1);
q2m = X(2);
q1mdot = X(3);
q2mdot = X(4);

if(t>=5 && t<=5.1)
    F = [1;1];
else
    F = [0;0];
end
    

firstm  = -inv(inertiam(q1m,q2m))*coriolism(q1m,q2m,q1mdot,q2mdot)*[q1mdot;q2mdot];
secondm =  inv(inertiam(q1m,q2m))*transpose(jacobm(q1m,q2m))*F;

finalm = firstm + secondm;

Dyn(1,1) = q1mdot;
Dyn(2,1) = q2mdot;
Dyn(3,1) = finalm(1,1);
Dyn(4,1) = finalm(2,1);

end

