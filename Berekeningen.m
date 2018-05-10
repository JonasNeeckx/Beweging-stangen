out = load('halfcycloiden.mat');
F_inert = out.normalforce_acc;
F_load = out.normalforce_load;
Pressure_Angle = out.pressure_angle;
rpm = out.rpm;
Mass = out.mass;
Omega_grad = rpm*360/60;
Omega_rad = rpm*2*pi/60;



