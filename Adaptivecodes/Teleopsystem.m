clear all;
%
%robot1
%{
massm = [3.14 2.26];
Im = [0.16 0.07];
lm = [1.04 0.96];
%}
%
massmcap = [3.00 2.00];
Imcap = [0.10 0.10];
lmcap = [0.8 1.2];
%}
%{
massmcap = [5.00 1.00];
Imcap = [0.60 0.80];
lmcap = [2 2];
%}
%{
massmcap = [1 1];
Imcap = [1 1];
lmcap = [1 1];
%}
%robot2
%{
masss = [2.54 1.82];
Is = [0.12 0.06];
ls = [1.28 0.84];
%}
%
massscap = [2.20 2.00];
Iscap = [0.2 0.1];
lscap = [1.1 0.5];
%}
%{
massscap = [4 3];
Iscap = [1 0.5];
lscap = [2 1];
%}
%{
massscap = [1 1];
Iscap = [1 1];
lscap = [1 1];
%}

thetadm1 =  massmcap(1)*(lmcap(1)/2)^2 + massmcap(2)*(lmcap(1)^2 + (lmcap(2)/2)^2 )+ Imcap(1)  + Imcap(2);
thetadm2 = massmcap(2)*lmcap(1)*lmcap(2)/2;
thetadm3 = massmcap(2)*(lmcap(2)/2)^2 + Imcap(2);

thetakm1 = lmcap(1);
thetakm2 = lmcap(2);


thetads1 =  massscap(1)*(lscap(1)/2)^2 + massscap(2)*(lscap(1)^2 + (lscap(2)/2)^2 )+ Iscap(1)  + Iscap(2);
thetads2 = massscap(2)*lscap(1)*lscap(2)/2;
thetads3 = massscap(2)*(lscap(2)/2)^2 + Iscap(2);

thetaks1 = lscap(1);
thetaks2 = lscap(2);



X0 = [ 0.8; 1.2; 0 ; 0; thetadm1; thetadm2;thetadm3;thetakm1; thetakm2; 1.0; 0.3; 0 ; 0; thetads1; thetads2;thetads3;thetaks1; thetaks2];
T = [0 60];
time = [0:1:T(2)];

[t,X] = ode15s(@(t,x) dynamics(t,x),T,X0);

%{
figure(1)
plot(t,X(:,1))

figure(2)
plot(t,X(:,5))
%legend('q1m','q1s')

figure(3)
plot(t,X(:,2))
%hold on

figure(4)
plot(t,X(:,6))
%legend('q2m','q2s')

figure(5)
plot(t, X(:,1) - X(:,5))

figure(6)
plot(t, X(:,2) - X(:,6))
%}

for i = 1:length(t)
    poss(i,:) = transs(X(i,10),X(i,11)) ;
    posm(i,:) = transm(X(i,1),X(i,2));
    
    vels(i,:) = (jacobs(X(i,10),X(i,11))*[X(i,12);X(i,13)]);
    velm(i,:) = (jacobm(X(i,1),X(i,2))*[X(i,3);X(i,4)]);
    
    Fh(i,:) = F_h(t(i),X(i,:));
    Fe(i,:) = F_e(t(i),X(i,:));
    
    thetadm(i,:) = [X(i,5);X(i,6);X(i,7)];
    thetads(i,:) = [X(i,14);X(i,15);X(i,16)];
    
    thetakm(i,:) = [X(i,8);X(i,9)];
    thetaks(i,:) = [X(i,17);X(i,18)];
    
    
    dyn(i,:) = dynamics(t(i),X(i,:));
    %Fconm(i,:) = inv(transpose(jacobm(X(i,1),X(i,2))))*(inertiam(X(i,1),X(i,2))*[dyn(i,3);dyn(i,4)] + coriolism(X(i,1),X(i,2),X(i,3),X(i,4))*[X(i,3);X(i,4)] - transpose(jacobm(X(i,1),X(i,2)))*F_h(t(i),X(i,:)));
    %Fcons(i,:) = inv(transpose(jacobs(X(i,5),X(i,6))))*(inertias(X(i,5),X(i,6))*[dyn(i,7);dyn(i,8)] + corioliss(X(i,5),X(i,6),X(i,7),X(i,8))*[X(i,7);X(i,8)] + transpose(jacobs(X(i,5),X(i,6)))*F_e(t(i),X(i,:)));
    
   Fconm(i,:) = inv(transpose(jacobm(X(i,1),X(i,2))))*(inertiam(X(i,1),X(i,2))*[dyn(i,3);dyn(i,4)] + inv(transpose(jacobm(X(i,1),X(i,2))))*coriolism(X(i,1),X(i,2),X(i,3),X(i,4))*[X(i,3);X(i,4)] - F_h(t(i),X(i,:)));
    Fcons(i,:) = inv(transpose(jacobs(X(i,10),X(i,11))))*(inertias(X(i,10),X(i,11))*[dyn(i,12);dyn(i,13)] + inv(transpose(jacobs(X(i,10),X(i,11))))*corioliss(X(i,10),X(i,11),X(i,12),X(i,13))*[X(i,12);X(i,13)] + F_e(t(i),X(i,:)));
    
end


figure(1)
plot(t,posm(:,1), 'b')
hold on
plot(t,poss(:,1),'b-.')
legend('Master Position','Slave Position')
xlabel('Time(s)')
ylabel('X Position (m)')


figure(2)
plot(t,posm(:,2),'b')
hold on
plot(t,poss(:,2),'b-.')
legend('Master Position','Slave Position')
xlabel('Time(s)')
ylabel('Y Position (m)')


figure(3)
plot(t, posm(:,1) - poss(:,1),'b')
%legend('Error in X')
xlabel('Time(s)')
ylabel('Error in Position (m)')
hold on
plot(t, posm(:,2) - poss(:,2),'b-.')
legend('Error in X','Error in Y')
xlabel('Time(s)')
%ylabel('Error in Y Position (m)')

figure(4)
plot(t, Fh(:,1),'b')
hold on
plot(t, Fe(:,1),'b-.')

ylabel('X axis Force (N)')
xlabel('Time(s)')
legend('Human Force Input','Environment Force')

figure(5)
plot(t, Fh(:,2),'b')
hold on
plot(t, Fe(:,2),'b-.')

ylabel('Y axis Force (N)')
xlabel('Time(s)')
legend('Human Force Input','Environment Force')

figure(6)
plot(t,thetadm(:,1), 'b')
hold on
plot(t,thetadm(:,2),'b-.')
plot(t,thetadm(:,3),'b--')
plot(time,4.0442*ones(61), 'g')
plot(time,1.1282*ones(61),'g-.')
plot(time,0.5907*ones(61),'g--')

legend('thetadm1','thetadm2','thetadm3', 'thetadm1star','thetadm2star','thetadm3star')
xlabel('Time(s)')
ylabel('Dynamic Parameters')

figure(7)
plot(t,thetads(:,1), 'r')
hold on
plot(t,thetads(:,2),'r-.')
plot(t,thetads(:,3),'r--')
plot(time,4.5233*ones(61), 'k')
plot(time,0.9784*ones(61),'k-.')
plot(time,0.3810*ones(61),'k--')

legend('thetads1','thetads2','thetads3','thetads1star','thetads2star','thetads3star')
xlabel('Time(s)')
ylabel('Dynamic Parameters')

figure(8)
plot(t,thetakm(:,1), 'b')
hold on
plot(t,thetakm(:,2),'b-.')

plot(time,1.04*ones(61), 'g')
plot(time,0.96*ones(61),'g-.')
xlabel('Time(s)')
ylabel('Kinematic Parameters')
legend('thetakm1','thetakm2','thetakm1star','thetakm2star')

figure(9)
plot(t,thetaks(:,1), 'r')
hold on
plot(t,thetaks(:,2),'r-.')
plot(time,1.28*ones(61), 'k')
plot(time,0.84*ones(61),'k-.')

legend('thetaks1','thetaks2','thetaks1star','thetaks2star')
xlabel('Time(s)')
ylabel('Kinematic Parameters')

%{
figure(10)

plot(t,Fh(:,2),'b')
hold on
plot(t,Fe(:,2),'b-.')
xlabel('Time(s)')
ylabel('Y axis Force (N)')
legend('Human Force Input','Environment Force')
%}

figure(11)
plot(t, Fconm(:,1),'b')
hold on
plot(t, Fcons(:,1),'b-.')

ylabel('X axis Force(N)')
xlabel('Time(s)')
legend('Master Control Force','Slave Control Force')

figure(12)
plot(t,Fconm(:,2),'b')
hold on
plot(t,Fcons(:,2),'b-.')
xlabel('Time(s)')
ylabel('Y axis Force(N)')
legend('Master Control Force','Slave Control Force')

%}