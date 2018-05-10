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
R0 = bcr + rof;

%% Geometry of the follower
%base circle and radius of the follower

%excentricity
e_optimal = optimalExcentricity(S,V,R0);