function [F] = F_h(t,X)
%FH Summary of this function goes here
%   Detailed explanation goes here
q1m = X(1);
q2m = X(2);
q1mdot = X(3);
q2mdot = X(4);
q1s = X(5);
q2s = X(6);
q1sdot = X(7);
q2sdot = X(8);

Rset1 = [0.9; 1.4];
Rset2 = [1.2; 1.1];
K =  75  ;
C =    50  ;

%K =  30  ;
%C =    15  ;


if(t<15)
    F = [0;0];
elseif (t<30)
    F = K*(Rset1-transm(q1m,q2m)) - C*(jacobm(q1m,q2m)*[q1mdot;q2mdot]);
elseif (t<45)
    F = K*(Rset2-transm(q1m,q2m)) - C*(jacobm(q1m,q2m)*[q1mdot;q2mdot]);
else
    F = [0;0];
end

end

