clear all;
%{
%robot1

massm = [3.14 2.26];
Im = [0.16 0.07];
lm = [1.04 0.96];


%robot2

masss = [2.54 1.82];
Is = [0.12 0.06];
ls = [1.28 0.84];
%}

X0 = [ 0.8; 1.2; 0 ; 0; 1.0; 0.3; 0 ; 0];
T = [0 60];

[t,X] = ode45(@(t,x) dynamics(t,x),T,X0);

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
    poss(i,:) = transs(X(i,5),X(i,6)) ;
    posm(i,:) = transm(X(i,1),X(i,2));
    
    vels(i,:) = (jacobs(X(i,5),X(i,6))*[X(i,7);X(i,8)]);
    velm(i,:) = (jacobm(X(i,1),X(i,2))*[X(i,3);X(i,4)]);
    
    Fh(i,:) = F_h(t(i),X(i,:));
    Fe(i,:) = F_e(t(i),X(i,:));
    
    
    dyn(i,:) = dynamics(t(i),X(i,:));
    %Fconm(i,:) = inv(transpose(jacobm(X(i,1),X(i,2))))*(inertiam(X(i,1),X(i,2))*[dyn(i,3);dyn(i,4)] + coriolism(X(i,1),X(i,2),X(i,3),X(i,4))*[X(i,3);X(i,4)] - transpose(jacobm(X(i,1),X(i,2)))*F_h(t(i),X(i,:)));
    %Fcons(i,:) = inv(transpose(jacobs(X(i,5),X(i,6))))*(inertias(X(i,5),X(i,6))*[dyn(i,7);dyn(i,8)] + corioliss(X(i,5),X(i,6),X(i,7),X(i,8))*[X(i,7);X(i,8)] + transpose(jacobs(X(i,5),X(i,6)))*F_e(t(i),X(i,:)));
    
    Fconm(i,:) = inv(transpose(jacobm(X(i,1),X(i,2))))*(inertiam(X(i,1),X(i,2))*[dyn(i,3);dyn(i,4)] + inv(transpose(jacobm(X(i,1),X(i,2))))*coriolism(X(i,1),X(i,2),X(i,3),X(i,4))*[X(i,3);X(i,4)] - F_h(t(i),X(i,:)));
    Fcons(i,:) = inv(transpose(jacobs(X(i,5),X(i,6))))*(inertias(X(i,5),X(i,6))*[dyn(i,7);dyn(i,8)] + inv(transpose(jacobs(X(i,5),X(i,6))))*corioliss(X(i,5),X(i,6),X(i,7),X(i,8))*[X(i,7);X(i,8)] + F_e(t(i),X(i,:)));
    
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
plot(t,Fh(:,2),'b')
hold on
plot(t,Fe(:,2),'b-.')
xlabel('Time(s)')
ylabel('Y axis Force (N)')
legend('Human Force Input','Environment Force')

figure(6)
plot(t, Fconm(:,1),'b')
hold on
plot(t, Fcons(:,1),'b-.')

ylabel('X axis Force(N)')
xlabel('Time(s)')
legend('Master Control Force','Slave Control Force')

figure(7)
plot(t,Fconm(:,2),'b')
hold on
plot(t,Fcons(:,2),'b-.')
xlabel('Time(s)')
ylabel('Y axis Force(N)')
legend('Master Control Force','Slave Control Force')

