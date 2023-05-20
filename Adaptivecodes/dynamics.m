function [Dyn] = dynamics(t,X) 
%DYNAMICS Summary of this function goes here
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

Fh = F_h(t,X);
Fe = F_e(t,X);


%
lambdam =  0.8   ;
kmm =   1    ;
krm =   4     ;
kdm =   1     ;
kpm = 0.5;
lambdas = 0.8    ;
kms =   1    ;
krs =   4     ;
kds =   1     ;
kps = 0.5;
taukm = 1*eye(2,2);
tauks = 1*eye(2,2);
taudm = 1*eye(3,3);
tauds = 1*eye(3,3);


%}

%{
lambdam =  0.5*0.8   ;
kmm =   0.5*1   ;
krm =   0.5*4     ;
kdm =   0.5*1     ;
lambdas = 0.5*0.8    ;
kms =   0.5*1    ;
krs =   0.5*4     ;
kds =    0.5*1    ;
taukm = 0.3*eye(2,2);
tauks = 0.3*eye(2,2);
taudm = 0.3*eye(3,3);
tauds = 0.3*eye(3,3);
%}


em = transs(q1s,q2s) - transm(q1m,q2m);
es = transm(q1m,q2m) - transs(q1s,q2s);

emdot = jacobs(q1s,q2s)*[q1sdot; q2sdot] - jacobm(q1m, q2m)*[q1mdot; q2mdot];
esdot = jacobm(q1m, q2m)*[q1mdot; q2mdot] - jacobs(q1s,q2s)*[q1sdot; q2sdot];

qrmdot = inv(jacobmcap(q1m,q2m,thetakcap1m,thetakcap2m))*lambdam*em;
qrsdot = inv(jacobmcap(q1s,q2s,thetakcap1s,thetakcap2s))*lambdas*es;

qrmddot = -inv(jacobmcap(q1m,q2m,thetakcap1m,thetakcap2m))*jacobmcapdot(q1m,q2m,q1mdot, q2mdot, thetakcap1m,thetakcap2m)*inv(jacobmcap(q1m,q2m,thetakcap1m,thetakcap2m))*lambdam*em + inv(jacobmcap(q1m,q2m,thetakcap1m,thetakcap2m))*lambdam*emdot;  
qrsddot = -inv(jacobmcap(q1s,q2s,thetakcap1s,thetakcap2s))*jacobmcapdot(q1s,q2s,q1sdot, q2sdot, thetakcap1s,thetakcap2s)*inv(jacobmcap(q1s,q2s,thetakcap1s,thetakcap2s))*lambdas*es + inv(jacobmcap(q1s,q2s,thetakcap1s,thetakcap2s))*lambdas*esdot;  


sm = [q1mdot;q2mdot] - qrmdot;
ss = [q1sdot;q2sdot] - qrsdot;

rm = inv(jacobmcap(q1m,q2m,thetakcap1m,thetakcap2m))*sm;
rs = inv(jacobmcap(q1s,q2s,thetakcap1s,thetakcap2s))*ss;

%{
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
%}

taum = Ydm(q1m,q2m,q1mdot,q2mdot,qrmdot(1,1),qrmdot(2,1),qrmddot(1,1),qrmddot(2,1))*[thetadcap1m;thetadcap2m;thetadcap3m] - kmm*sm - krm*transpose(jacobmcap(q1m,q2m,thetakcap1m,thetakcap2m))*rm + kdm*transpose(jacobmcap(q1m,q2m,thetakcap1m,thetakcap2m))*emdot;
taus = Ydm(q1s,q2s,q1sdot,q2sdot,qrsdot(1,1),qrsdot(2,1),qrsddot(1,1),qrsddot(2,1))*[thetadcap1s;thetadcap2s;thetadcap3s] - kms*ss - krs*transpose(jacobmcap(q1s,q2s,thetakcap1s,thetakcap2s))*rs + kds*transpose(jacobmcap(q1s,q2s,thetakcap1s,thetakcap2s))*esdot;

qmddot = inv(inertiam(q1m,q2m))*(taum - coriolism(q1m,q2m,q1mdot,q2mdot)*[q1mdot;q2mdot] + transpose(jacobm(q1m,q2m))*Fh) ;
qsddot = inv(inertias(q1s,q2s))*(taus - corioliss(q1s,q2s,q1sdot,q2sdot)*[q1sdot;q2sdot] - transpose(jacobs(q1s,q2s))*Fe) ;

thetadmdot = -taudm*(transpose(Ydm(q1m,q2m,q1mdot,q2mdot,qrmdot(1,1),qrmdot(2,1),qrmddot(1,1),qrmddot(2,1)))*sm + 0.01*[thetadcap1m;thetadcap2m;thetadcap3m]);
thetadsdot = -tauds*(transpose(Ydm(q1s,q2s,q1sdot,q2sdot,qrsdot(1,1),qrsdot(2,1),qrsddot(1,1),qrsddot(2,1)))*ss + 0.01*[thetadcap1s;thetadcap2s;thetadcap3s]);

deltaXm = transmcap(q1m,q2m,thetakcap1m,thetakcap2m) - transm(q1m,q2m);
deltaXs = transmcap(q1s,q2s,thetakcap1s,thetakcap2s) - transs(q1s,q2s);

thetakmdot = -taukm*(kpm*transpose(Ykm(q1m,q2m,q1mdot,q2mdot))*deltaXm + kdm*transpose(Ykm(q1m,q2m,q1mdot,q2mdot))*emdot + 0.01*[thetakcap1m;thetakcap2m]);
thetaksdot = -tauks*(kps*transpose(Ykm(q1s,q2s,q1sdot,q2sdot))*deltaXs + kds*transpose(Ykm(q1s,q2s,q1sdot,q2sdot))*esdot + 0.01*[thetakcap1s;thetakcap2s]);

Dyn(1,1) = q1mdot;
Dyn(2,1) = q2mdot;
Dyn(3,1) = qmddot(1,1) ;
Dyn(4,1) = qmddot(2,1);
Dyn(5,1) = thetadmdot(1,1);
Dyn(6,1) = thetadmdot(2,1);
Dyn(7,1) = thetadmdot(3,1);
Dyn(8,1) = thetakmdot(1,1);
Dyn(9,1) = thetakmdot(2,1);


Dyn(10,1) = q1sdot;
Dyn(11,1) = q2sdot;
Dyn(12,1) = qsddot(1,1);
Dyn(13,1) = qsddot(2,1);
Dyn(14,1) = thetadsdot(1,1);
Dyn(15,1) = thetadsdot(2,1);
Dyn(16,1) = thetadsdot(3,1);
Dyn(17,1) = thetaksdot(1,1);
Dyn(18,1) = thetaksdot(2,1);

end

