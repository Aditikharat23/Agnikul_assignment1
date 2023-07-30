clear all
close all
m1 = 10e15;
r1  = [0,0];
v1 = [10,0];

m2 = 20;
r2 = [0,10];
v2 = [-10,0];

IC = [r1,r2,v1,v2];
tend = 500;
[t,state] = RK4.conststep(IC,0.5,@tbp,tend);
%% 


[t2,state2] = ode45(@tbp,[0 tend],IC);
%% Plotting
figure(1)
plot(state(:,1),state(:,2))
hold on
plot(state(:,3),state(:,4))
legend('m1','m2')
title('Constant Step size of 0.5s')
figure(2)
plot(state2(:,1),state2(:,2))
hold on
plot(state2(:,3),state2(:,4))
legend('m1','m2')
title('Variable step size ode45')
%% 

function dstatedt = tbp(t,state)
%State = [r1,r2,v1,v2];
m1 =10e23;
m2=20e23;
G = 6.67430*10^-11;

R = state(3:4) - state(1:2);
a1 = G*m2*R/norm(R)^3;
a2 = -G*m1*R/norm(R)^3;

dstatedt = [state(5:8)',a1(1),a1(2),a2(1),a2(2)]';
end


