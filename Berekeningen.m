%% Load parameters
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

k = 331.8*cosd(25.33)/14.99;
[Springconstant_optimal, Fv0_optimal,optimal_N] = spring(S, F_load, F_inert, Pressure_Angle)


%% Geometry of the follower
%base circle and radius of the follower

%excentricity
e_optimal = optimalExcentricity(S,V,R0);