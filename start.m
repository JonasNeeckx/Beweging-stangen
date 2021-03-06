%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Analysis of a 9 bar linkage system
% 
% Jonas Neeckx
% Nicolas Heintz
% 
% Based on the work of:
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_9bar = 1;           % draw figures of kinematic analysis if 1
fig_kin_check = 1;          % draw figures of kinematic checks if 1
fig_dyn_9bar = 0;           % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
r2 = 0.00375;
r3 = 0.00786;
r5 = 0.00800;
r6 = 0.00450;
r7 = 0.00800;
r9 = 0.00800;
r10 = 0.00450;
r11 = 0.00800;
r14y = 0.00691;
r14x = 0.00350;
r47y = 0.00350;
r18x = 0.00350;
r18y = 0.00691;
r811y = 0.0035;

% dynamic parameters, defined in a local frame on each of the bars.
% the centre of gravity of bar 8 and 4 are calculated with respect to their
% maximal lenght, instead of the varying length r8 and r4.

rmax4 = sqrt(r14x^2 + (r14y - r2 - r3)^2);
rmax8 = sqrt(r18x^2 + (r18y - r2 - r3)^2);
rwing = 0.06;

X2 = r2/2;               % X coordinates of cog (centre of gravity)
X3 = r3/2;
X4 = rmax4/2;
X5 = r5/2;
X6 = r6/2;
X7 = r7/2;
X8 = rmax8/2;
X9 = r9/2;
X10 = r10/2;
X11 = r11/2;
Xwing = 0;

Y2 = 0;                  % Y coordinates of cog
Y3 = 0;
Y4 = 0;
Y5 = 0;
Y6 = -rwing^2/(2*(rwing+r6));
Y7 = 0;
Y8 = 0;
Y9 = 0;
Y10 = rwing^2/(2*(rwing+r10));
Y11 = 0;
Ywing = rwing/2;

dbar = 2*10^-3;
rho_link = 1800*pi*dbar^2/4;
m2 = r2*rho_link;
m3 = r3*rho_link;
m4 = rmax4*rho_link;
m5 = r5*rho_link;
m6 = r6*rho_link;
m7 = r7*rho_link;
m8 = rmax8*rho_link;
m9 = r9*rho_link;
m10 = r10*rho_link;
m11 = r11*rho_link;
mwing = rwing*rho_link;

%J4, J5, J8 and J9 are defined with respect to their fixed point
%For J6 and J10 the inertia of the wing is taken into account
J2 = m2*r2^2/12;
J3 = m3*r3^2/12;
J4 = m4*rmax4^2/12 + m4*X4^2;
J5 = m5*r5^2/12 + m5*X5^2;
J6 = m6*r6^2/12 + m6*Y6^2 + mwing*rwing^2/12 + mwing*(-rwing/2-Y6)^2;
J7 = m7*r7^2/12;
J8 = m8*rmax8^2/12 + m8*X8^2;
J9 = m9*r9^2/12 + m9*X9^2;
J10 = m10*r10^2/12 + m10*Y10^2 + mwing*rwing^2/12 + mwing*(rwing/2-Y10)^2;
J11 = m11*r11^2/12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis
phi3_init = 60*pi/180;
phi4_init = 0;    % initial condition for first step of position analysis with fsolve
phi6_init = 90*pi/180;  % VERY IMPORTANT because it determines which branch of the mechanism you're in
phi7_init = 210*pi/180;
phi8_init = pi;
phi10_init = 90*pi/180;
phi11_init = 330*pi/180;
r13_init = 0.00691;
r4_init = 0.00350;
r8_init = 0.00350; 

t_begin = 0;                   % start time of simulation
t_end = 4/37;                    % end time of simulation
Ts = 0.0001;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 37*2*pi;                   % omega = 145
phi2=omega*t+(0*pi/180);
tijdsvec = size(t);
dphi2=ones(tijdsvec(1),1).*omega; % Always omega
ddphi2=zeros(tijdsvec(1),1); % Always zero

% calculation of the kinematics (see kin_4bar.m)
[phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,dphi3,dphi4,dphi6,dphi7,dphi8,dphi10,dphi11,ddphi3,ddphi4,ddphi6,ddphi7,ddphi8,ddphi10,ddphi11,r8,dr8,ddr8,r13,dr13,ddr13,r4,dr4,ddr4] = kinematics_9bar(r2,r3,r5,r6,r7,r9,r10,r11,r14x,r14y,r47y,r18x,r18y,r811y,phi2,dphi2,ddphi2,phi3_init,phi4_init,phi6_init,phi7_init,phi8_init,phi10_init,phi11_init,r13_init,r4_init,r8_init,t,fig_kin_9bar,fig_kin_check);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
[F_A_x, F_A_y, F_23_x, F_23_y, F_C_x, F_34, F_38, F_D_x, F_D_y, F_56_x, F_56_y, F_67_x, F_67_y,...
    F_G_x, F_G_y, F_H_x, F_H_y, F_910_x, F_910_y, F_1011_x, F_1011_y, F_K_x, F_K_y, M_A] = ...
    dynamics_9bar(phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,...
dphi2,dphi3,dphi4, dphi6, dphi7, dphi8, dphi10, dphi11,...
ddphi2,ddphi3,ddphi4,ddphi6, ddphi7, ddphi8, ddphi10, ddphi11,...
r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,rwing,rmax4,rmax8,...
m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,mwing,...
X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,Xwing,...
Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Ywing,...
J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,t,fig_dyn_9bar);

[F_A_x_bal, F_A_y_bal, F_23_x_bal, F_23_y_bal, F_C_x_bal, F_34_bal, F_38_bal, F_D_x_bal, F_D_y_bal, F_56_x_bal, F_56_y_bal, F_67_x_bal, F_67_y_bal,...
    F_G_x_bal, F_G_y_bal, F_H_x_bal, F_H_y_bal, F_910_x_bal, F_910_y_bal, F_1011_x_bal, F_1011_y_bal, F_K_x_bal, F_K_y_bal, M_A_bal] = ...
    dynamics_9bar_balanced(phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,...
dphi2,dphi3,dphi4, dphi6, dphi7, dphi8, dphi10, dphi11,...
ddphi2,ddphi3,ddphi4,ddphi6, ddphi7, ddphi8, ddphi10, ddphi11,...
r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,rwing,rmax4,rmax8,...
m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,mwing,...
X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,Xwing,...
Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Ywing,...
J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,t,fig_dyn_9bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)

