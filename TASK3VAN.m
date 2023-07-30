%%Initial Condtition
clear all
close all
ic = [0.02 0];
tend = 30;
m   = 0.3;
             b  = 0.77;
             k=0.11;
             root  = roots([m b k]);

% syms C1 C2 
% eqn1 = C1 + C2 == ic(1);
% eqn2 = root(1)*C1 + root(2)*C2 == -ic(2);
% 
% [A,B] = equationsToMatrix([eqn1, eqn2], [C1, C2]);
% X = linsolve(A,B)
%% Constant step size

h1 = 1; %%Step size of 1s
h2 = 0.5; %Step size of 0.5s
h3 = 0.2; %Step size of 0.1s

[t1,state1] = R_k.conststep(ic,h1,@R_k.dmp,tend);
[t2,state2] = R_k.conststep(ic,h2,@R_k.dmp,tend);
[t3,state3] = R_k.conststep(ic,h3,@R_k.dmp,tend);

%% Analtical sol
for i =1:length(t1)
ref1(i,:) = R_k.sol_dmp(t1(i));
error1(i,:) = abs(ref1(i,:)-state1(i,1:2)); %Error in x and y 
end
for i =1:length(t2)
ref2(i,:) = R_k.sol_dmp(t2(i));
error2(i,:) =abs(ref2(i,:)-state2(i,1:2)); %Error in x and y 
end
for i =1:length(t3)
ref3(i,:) = R_k.sol_dmp(t3(i));
error3(i,:) = abs(ref3(i,:)-state3(i,1:2)); %Error in x and y 
end

%% %% Adaptive step size
[t4,state4,h_track] = R_k.stephalving45(ic,10^-5,@R_k.dmp,tend,1);
for i =1:length(t4)
ref4(i,:) = R_k.sol_dmp(t4(i));
error4(i,:) = abs(ref4(i,:)-state4(i,1:2)); %Error in x and y 
end

%% PLOTTING
figure(1)
plot(t1,error1(:,1))
hold on 
plot(t2,error2(:,1))
hold on 
plot(t3,error3(:,1))
hold on 
plot(t4,error4(:,1))

title('Error in position')
legend('h=1','h=0.5','h=0.2','h =variable')
xlabel('Time (s)')
ylabel('Error')

figure(2)
plot(t4,state4(:,1))
hold on 
plot(t4,ref4(:,1))
xlabel('Time(s)')
ylabel('X (m)')
title('Mass Damper System motion')




