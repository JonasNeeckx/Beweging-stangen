%% Load parameters
clear; clc; 
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
R0 = bcr + rof;Omega_grad = rpm*360/60;
Omega_rad = rpm*2*pi/60;
normalForce = out.normalforce_tot;
xcam = out.xcam;
ycam = out.ycam;

k = 331.8*cosd(25.33)/14.99;

out_no_e = load('halfcycloiden_no_exc');
Pressure_Angle_no_e = out_no_e.pressure_angle;
V_no_e = out_no_e.V;
normalForce_no_e = out_no_e.normalforce_tot;

%% Geometry of the follower
%base circle and radius of the follower

%excentricity
e_optimal = optimalExcentricity(S,V,R0);

%% Rigid-body forces
%dimensionalisation of the spring

%instantaneous power
instantanious_power = instantaniousPower(normalForce,Pressure_Angle,rpm,V);
instantanious_power_no_e = instantaniousPower(normalForce_no_e,Pressure_Angle_no_e,rpm,V_no_e);

%average power
average_power = mean(instantanious_power);



