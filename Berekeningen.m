%% Load parameters
clear; clc; close all;
out = load('halfcycloiden.mat');
F_inert = out.normalforce_acc;
F_load = out.normalforce_load;
Pressure_Angle = out.pressure_angle;
rpm = out.rpm;
Mass = out.mass;
S = out.S;
V = out.V;
bcr = out.bcr;
rof = out.rof;
R0 = bcr + rof;
Omega_grad = rpm*360/60;
Omega_rad = rpm*2*pi/60;
normalForce = out.normalforce_tot;
xcam = out.xcam;
ycam = out.ycam;
k = 331.8*cosd(25.33)/14.99;

out_no_e = load('halfcycloiden_no_exc');
Pressure_Angle_no_e = out_no_e.pressure_angle;
V_no_e = out_no_e.V;
normalForce_no_e = out_no_e.normalforce_tot;
zeta = 0.1;
theta1 = out.theta;
theta = theta1;

%% Geometry of the follower
%base circle and radius of the follower

%excentricity
e_optimal = optimalExcentricity(S,V,R0);

%% Rigid-body forces
%dimensionalisation of the spring
[ks, Fv0_optimal] = spring(S, F_load, F_inert, Pressure_Angle);

%instantaneous power
instantaneous_power = instantaniousPower(normalForce,Pressure_Angle,rpm,V);
instantaneous_power_no_e = instantaniousPower(normalForce_no_e,Pressure_Angle_no_e,rpm,V_no_e);

%average power
average_power = mean(instantaneous_power);

plot_power(average_power, instantaneous_power, instantaneous_power_no_e, theta)

%instantanious torque
instantanious_torque = instantaneous_power./Omega_rad;

%average torque
average_torque = mean(instantanious_torque);

% I flywheel
[I_flywheel,mean_index] = flywheel(instantanious_torque,average_torque,Omega_rad);

% Speed variation
w = speed_variation(average_torque, instantanious_torque,I_flywheel, mean_index, Omega_rad);

%% Dynamics of a deformable follower
%single rise
epsilon = Single_Rise(zeta,ks,theta1,Omega_rad,Mass,S,normalForce);

