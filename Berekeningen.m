out = load('halfcycloiden.mat');
F_inert = out.normalforce_acc;
F_load = out.normalforce_load;
Pressure_Angle = out.pressure_angle;
rpm = out.rpm;
Mass = out.mass;
Omega_grad = rpm*360/60;
Omega_rad = rpm*2*pi/60;
Lift = out.S;
normalForce = out.normalforce_tot;

k = 331.8*cosd(25.33)/14.99;
[Springconstant_optimal, Fv0_optimal] = spring(Lift, F_load, F_inert, Pressure_Angle)


