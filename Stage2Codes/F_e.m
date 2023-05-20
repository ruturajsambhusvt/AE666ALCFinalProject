function [F] = F_e(t,X)
%FE Summary of this function goes here
%   Detailed explanation goes here

q1m = X(1);
q2m = X(2);
q1mdot = X(3);
q2mdot = X(4);
thetadcap1m = X(5);
thetadcap2m = X(6);
thetadcap3m = X(7);
thetakcap1m = X(8);
thetakcap2m = X(9);

q1s = X(10);
q2s = X(11);
q1sdot = X(12);
q2sdot = X(13);
thetadcap1s = X(14);
thetadcap2s = X(15);
thetadcap3s = X(16);
thetakcap1s = X(17);
thetakcap2s = X(18);

%Ke =  500    ;
%Ce =   0.1    ;
Ke =  80    ;
Ce =   0.1    ;
Xwall = 1;

positions = transs(q1s,q2s);
velocitys = jacobs(q1s,q2s)*[q1sdot;q2sdot];
if positions(1)>1
    F(1,1) = Ke*(positions(1)-Xwall) + Ce*velocitys(1);
else
    F(1,1) = 0;
    
end

F(2,1) = 0;
end

