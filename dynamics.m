function [Dyn] = dynamics(t,X) 
%DYNAMICS Summary of this function goes here
%   Detailed explanation goes here

q1m = X(1);
q2m = X(2);
q1mdot = X(3);
q2mdot = X(4);
q1s = X(5);
q2s = X(6);
q1sdot = X(7);
q2sdot = X(8);

Fh = F_h(t,X);
Fe = F_e(t,X);


%{
lambdam =  0.8   ;
kmm =   1    ;
krm =   4     ;
kdm =   1     ;
lambdas = 0.8    ;
kms =   1    ;
krs =   4     ;
kds =   1     ;
%}

%
lambdam =  0.5*0.8   ;
kmm =   0.5*1   ;
krm =   0.5*4     ;
kdm =   0.5*1     ;
lambdas = 0.5*0.8    ;
kms =   0.5*1    ;
krs =   0.5*4     ;
kds =    0.5*1    ;
%}

firstm = lambdam*(-inv(jacobm(q1m,q2m))*jacobmdot(q1m,q2m,q1mdot,q2mdot)*inv(jacobm(q1m,q2m)))*(transs(q1s,q2s)-transm(q1m,q2m));
secondm = lambdam*(inv(inertiam(q1m,q2m))*coriolism(q1m,q2m,q1mdot,q2mdot) *inv(jacobm(q1m,q2m)))*(transs(q1s,q2s)-transm(q1m,q2m));
thirdm = lambdam*kmm*(inv(inertiam(q1m,q2m))*inv(jacobm(q1m,q2m)))*(transs(q1s,q2s)-transm(q1m,q2m));
fourthm = krm*(inv(inertiam(q1m,q2m))*transpose(jacobm(q1m,q2m)))*(transs(q1s,q2s)-transm(q1m,q2m));

fifthm = lambdam*inv(jacobm(q1m,q2m))*((jacobs(q1s,q2s)*[q1sdot;q2sdot])- (jacobm(q1m,q2m)*[q1mdot;q2mdot]));
sixthm = -kdm*(inv(inertiam(q1m,q2m))*transpose(jacobm(q1m,q2m)))*((jacobs(q1s,q2s)*[q1sdot;q2sdot])- (jacobm(q1m,q2m)*[q1mdot;q2mdot]));

seventhm = -inv(inertiam(q1m,q2m))*coriolism(q1m,q2m,q1mdot,q2mdot)*[q1mdot;q2mdot];
eighthm = -kmm*inv(inertiam(q1m,q2m))*[q1mdot;q2mdot];
ninthm = -krm*(inv(inertiam(q1m,q2m))*transpose(jacobm(q1m,q2m))*jacobm(q1m,q2m))*[q1mdot;q2mdot];
tenthm = transpose(jacobm(q1m,q2m))*Fh;


finalm = firstm + secondm + thirdm + fourthm + fifthm + sixthm + seventhm + eighthm + ninthm + tenthm; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firsts = lambdas*(-inv(jacobs(q1s,q2s))*jacobsdot(q1s,q2s,q1sdot,q2sdot)*inv(jacobs(q1s,q2s)))*(transm(q1m,q2m)-transs(q1s,q2s));
seconds = lambdas*(inv(inertias(q1s,q2s))*corioliss(q1s,q2s,q1sdot,q2sdot) *inv(jacobs(q1s,q2s)))*(transm(q1m,q2m)-transs(q1s,q2s));
thirds = lambdas*kms*(inv(inertias(q1s,q2s))*inv(jacobs(q1s,q2s)))*(transm(q1m,q2m)-transs(q1s,q2s));
fourths = krs*(inv(inertias(q1s,q2s))*transpose(jacobs(q1s,q2s)))*(transm(q1m,q2m)-transs(q1s,q2s));

fifths = lambdas*inv(jacobs(q1s,q2s))*((jacobm(q1m,q2m)*[q1mdot;q2mdot])- (jacobs(q1s,q2s)*[q1sdot;q2sdot]));
sixths = -kds*(inv(inertias(q1s,q2s))*transpose(jacobs(q1s,q2s)))*((jacobm(q1m,q2m)*[q1mdot;q2mdot])- (jacobs(q1s,q2s)*[q1sdot;q2sdot]));

sevenths = -inv(inertias(q1s,q2s))*corioliss(q1s,q2s,q1sdot,q2sdot)*[q1sdot;q2sdot];
eighths = -kms*inv(inertias(q1s,q2s))*[q1sdot;q2sdot];
ninths = -krs*(inv(inertias(q1s,q2s))*transpose(jacobs(q1s,q2s))*jacobm(q1s,q2s))*[q1sdot;q2sdot];
tenths = -transpose(jacobs(q1s,q2s))*Fe;

finals = firsts + seconds + thirds + fourths + fifths + sixths + sevenths + eighths + ninths+ tenths ; 

Dyn(1,1) = q1mdot;
Dyn(2,1) = q2mdot;
Dyn(3,1) = finalm(1,1);
Dyn(4,1) = finalm(2,1);

Dyn(5,1) = q1sdot;
Dyn(6,1) = q2sdot;
Dyn(7,1) = finals(1,1);
Dyn(8,1) = finals(2,1);


end

