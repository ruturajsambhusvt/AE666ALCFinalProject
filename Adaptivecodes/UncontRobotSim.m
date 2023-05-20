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

X0 = [ 0.8; 1.2; 0 ; 0];
T = [0 100];

[t,X] = ode23(@(t,x) UnconRobot(t,x),T,X0);

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
    %poss(i,:) = transs(X(i,5),X(i,6)) ;
    posm(i,:) = transm(X(i,1),X(i,2));
    
    %vels(i,:) = (jacobs(X(i,5),X(i,6))*[X(i,7);X(i,8)]);
    velm(i,:) = (jacobm(X(i,1),X(i,2))*[X(i,3);X(i,4)]);
    
end


figure(1)
plot(t,posm(:,1),'k')
hold on
%plot(t,poss(:,1))
%legend('Xm')
xlabel('time(s)')
ylabel('Positon (m)')

%figure(2)
plot(t,posm(:,2),'r')
hold on
%plot(t,poss(:,2))
legend('Xm', 'Ym')
%xlabel('time(s)')
%ylabel('Y position (m)')

%{
figure(5)
plot(t, posm(:,1) - poss(:,1))
legend('Error in X')

figure(6)
plot(t, posm(:,2) - poss(:,2))
legend('Error in Y')
%}
