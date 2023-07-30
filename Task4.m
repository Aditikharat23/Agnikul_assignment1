%%Rocket Ascent Trajectory Code
%Constants
clear all
close all
global g0 Re h_grav_turn A Cd T Isp m_p rho0 m_s
g0 =9.80665; %m/s^2 Gravitational acceleration at sea level 
Re = 6378000; %Radius of earth
rho0 = 1.225; %Density at sea level
%% Given Parameters
T = 100000;             %Thrust [Newtons]
t2w = 1.4;              %Thrust to weight ratio
m_total = T/(t2w*g0);
h_grav_turn = 100; % Height for gravity turn [m]

%% Assumed parameters [Skylon]
d = 4; %Diameter [m]
Isp = 350;% Specific Impulse [s]
s = 0.02; %Structural Ratio
Cd = 0.5 ; %Drag Coefficient assumed to be constant

%% Calculating additional parameters based on given data and assumptions
m_p = m_total*(1-s); %Mass of propellant [kg]
m_dot = T/(Isp*g0); % Mass flow rate [kg/s]
t_burn = m_p/m_dot; %Approximate burn time assuming constant thrust[s]
m_s = s*m_total;
A = (pi*d^2)/4; %Cross sectional area[m^2]
%% Initial Conditions
v0 = 0; %Initial Velocity
x0 = 0; %Initial Downrange distance
h0 = 0; %Initial Altitude
fpa0 = deg2rad(89.999); %Flight path angle at gravity turn
m_0 = m_total; %Initial mass of the rocket


IC = [v0,x0,h0,fpa0,m_0];
tend = t_burn+15; %To account for throttle down during maxQ
%% Integrating launch ODE using inbuilt ode45
[t,state] = ode45(@ascent,[0 tend],IC);
%% Defining the state variables
v = state(:,1)/1000; %Velocity in km/s
x = state(:,2)/1000; %Downrange distance in km
h = state(:,3)/1000; %Height in km
fpa = rad2deg(state(:,4)); % Flight path angle in deg
m = state(:,5); %Mass remaining in kg

%
for i = 1:length(t)

[dum,rho] = atmosnrlmsise00(h(i)*1000,13.7259,80.2266,2023,1,0);
Rho(i) = rho(6);

 q(i) = 1/2*Rho(i)*(v(i)*1000)^2; %Dynamic pressure (Pa)

end

q = q/101325; %Pascal to atm
[maxQ,imax] = max(q);
tQ = t(imax);
vQ = v(imax); 
hQ = h(imax); 

%% 
figure(1)
subplot(2,2,1)
plot(t,h)
xlabel("Time(s)")
ylabel("Altitude (km)")
subplot(2,2,2)
plot(t,fpa)
xlabel("Time(s)")
ylabel("Flight Path Angle(Degrees)")
subplot(2,2,3)
plot(t,m)
xlabel("Time(s)")
ylabel("Mass (kg)")
subplot(2,2,4)
plot(t,v)
xlabel("Time(s)")
ylabel("Velocity (km/s)")

figure(2)
plot(h,q)
title('dynamic Pressure Variation with height')
xlabel('Height (km)')
ylabel('Dyanmic Pressure (atm)')
%% 
function dstatedt = ascent(t,state)
global g0 Re h_grav_turn A Cd T Isp m_p rho0 m_s
persistent flag

v = state(1);
x = state(2);
h = state(3);
fpa = state(4);
m = state(5);
%% Drag 
[dummy,Rho] = atmosnrlmsise00(h,45,-50,2023,4,0);rho = Rho(6);
% rho = rho0*exp(-h/7500);
% [dummy, dummy, dummy, rho] = atmosisa(h);
D = 0.5*rho*(v^2)*A*Cd;
%% Gravity
g = g0/(1 + h/Re)^2 ;%Variation of gravity with height
%% 
if m>m_s
    T = 100000;
else T=0;
end
if h < h_grav_turn

v_dot = T/m - D/m - g;
x_dot = 0;
h_dot = v;
fpa_dot = 0;

else

v_dot = T/m - D/m - g*sin(fpa);
fpa_dot = -(1/v)*(g - (v^2)/(Re + h))*cos(fpa)

x_dot = (Re/(Re + h))*v*cos(fpa); 
h_dot = v*sin(fpa);
end
m_dot = -T/(Isp*g0);

dstatedt = [v_dot,x_dot,h_dot,fpa_dot,m_dot]';
end






