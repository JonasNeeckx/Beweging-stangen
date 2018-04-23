%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
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
fig_kin_4bar = 1;        % draw figures of kinematic analysis if 1
fig_dyn_4bar = 0;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
r2 = 0.00375;
r3 = 0.00786;
r4 = 0.0035;
r5 = 0.00800;
r6 = 0.0045;
r7 = 0.00800;
r9 = 0.008;
r10 = 0.0045;
r11 = 0.008;
r14y = 0.00691;
r14x = 0.0035;
r47y = 0.0045;
r18x = 0.0035;
r18y = 0.00691;
r811y = 0.0045;

% dynamic parameters, defined in a local frame on each of the bars.
X2 = r2/2;               % X coordinates of cog (centre of gravity)
X3 = r3/3;
X4 = r4/2;
X5 = r5/2;
X6 = r6/2;
X7 = r7/2;
X8 = r8/2;
X9 = r9/2;
X10 = r10/2;
X11 = r11/2;

Y2 = 0;                  % Y coordinates of cog
Y3 = 0;
Y4 = 0;
Y5 = 0;
Y6 = 0;
Y7 = 0;
Y8 = 0;
Y9 = 0;
Y10 = 0;
Y11 = 0;

m2 = r2*0.5;
m3 = r3*0.5;
m4 = r4*0.5;
m5 = r5*0.5;
m6 = r6*0.5;
m7 = r7*0.5;
m9 = r9*0.5;
m10 = r10*0.5;
m11 = r11*0.5;

J2 = m2*r2^2/12;
J3 = m3*r3^2/12;
J4 = m4*r4^2/12;
J5 = m5*r5^2/12;
J6 = m6*r6^2/12;
J7 = m7*r7^2/12;
J9 = m9*r9^2/12;
J10 = m9*r9^2/12;
J11 = m9*r9^2/12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis
phi2_init = 5.8;
phi3_init = 1.6;    % initial condition for first step of position analysis with fsolve (phi3 and phi4)
phi5_init = 1.6;  % VERY IMPORTANT because it determines which branch of the mechanism you're in
phi6_init = 2.1;
phi8_init = 4.7;
phi10_init = 1.6;
phi11_init = 3.7;
r8_init = 0.035; 

t_begin = 0;                   % start time of simulation
t_end = 10;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 0.5;                   % omega = 145
phi1=omega*t+(60*pi/180);
tijdsvec = size(t);
dphi1=ones(tijdsvec(1),1).*omega; % ten allen tijde omega
ddphi1=zeros(tijdsvec(1),1); % ten allen tijde nul

% calculation of the kinematics (see kin_4bar.m)
[phi2,phi3,phi5,phi6,phi8,phi10,phi11,dphi2,dphi3,dphi5,dphi6,dphi8,dphi10,dphi11,ddphi2,ddphi3,ddphi5,ddphi6,ddphi8,ddphi10,ddphi11,r8,dr8,ddr8] = kinematics_4bar(r2,r3,r4,r5,r6,r7,r9,r10,r11,r14x,r14y,r47y,r18x,r18y,r811y,phi1,dphi1,ddphi1,phi2_init,phi3_init,phi5_init,phi6_init,phi8_init,phi10_init,phi11_init,r8_init,t,fig_kin_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
[F_P_x,F_Q_x,F_R_x,F_S_x,F_P_y,F_Q_y,F_R_y,F_S_y,M_P] = dynamics_4bar(phi2,phi3,phi4,phi5,phi6,phi8,phi9,phi10,phi11,dphi2,dphi3,dphi4,dphi5,dphi6,dphi8,dphi9,dphi10,dphi11,ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi8,ddphi9,ddphi10,ddphi11,r2,r3,r4, ...
  r5,r6,r7,r9,r10,r11,r14x,r14y,r18x,r18y,r47y,r811y,m2,m3,m4,m4,m5,m6,m7,m8,m9,m10,m11,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)

