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


function [F_A_x, F_A_y, F_23_x, F_23_y, F_C_x, F_34, F_38, F_D_x, F_D_y, F_56_x, F_56_y, F_67_x, F_67_y,...
    F_G_x, F_G_y, F_H_x, F_H_y, F_910_x, F_910_y, F_1011_x, F_1011_y, F_K_x, F_K_y, M_A] = ...
dynamics_4bar(phi1,phi2,phi3,phi4,phi5,phi6,phi8,phi9,phi10,phi11,...
dphi1,dphi2,dphi3, dphi5, dphi6, dphi8,  dphi10, dphi11,...
ddphi1,ddphi2,ddphi3,ddphi5, ddphi6, ddphi8,  ddphi10, ddphi11,...
r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,rwing, ...
m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,mwing, ...
X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,Xwing, ...
Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Ywing, ...
J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,t,fig_dyn_4bar)

% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.

%As angle phi4 and phi9 are equal to respectivily phi3 and phi8 apart from
%a constant value, their dphi and ddphi is not given.
dphi4 = dphi3;
ddphi4 = ddphi3;
dphi9 = dphi8;
ddphi9 = ddphi8;

%The masses and moments of inertia of bonded links are added
%together
J45 = J4 + J5;
J89 = J8 + J9;
m45 = m4 + m5;
m89 = m8 + m9;
m6wing = m6 + mwing;
m10wing = m10 + mwing;

%Thrust on the wing (simplified to a force in the middle of the wing)
F_wing = 0.12;

% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P
cog2_A_x = X2*cos(phi1);
cog2_A_y = X2*sin(phi1);
cog2_B_x = -X2*cos(phi1);
cog2_B_y = -X2*sin(phi1);
cog3_B_x = -X3*cos(phi2);
cog3_B_y = -X3*sin(phi2);
cog3_C_x = X3*cos(phi2);
cog3_C_y = X3*sin(phi2); 
cog6_E_x = cos(phi5)*X6 - sin(phi5)*Y6;
cog6_E_y = sin(phi5)*X6 + cos(phi5)*Y6;
cog6_F_x = -cos(phi5)*X6 - sin(phi5)*Y6;
cog6_F_y = -sin(phi5)*X6 + cos(phi5)*Y6;
cog7_F_x = cos(phi6)*X7;
cog7_F_y = sin(phi6)*X7;
cog7_G_x = -cos(phi6)*X7;
cog7_G_y = -sin(phi6)*X7;
cog10_I_x = cos(phi10)*X10 - sin(phi10)*Y10;
cog10_I_y = sin(phi10)*X10 + cos(phi10)*Y10;
cog10_J_x = -cos(phi10)*X10 - sin(phi10)*Y10;
cog10_J_y =	-sin(phi10)*X10 + cos(phi10)*Y10;
cog11_J_x = cos(phi11)*X11;
cog11_J_y = sin(phi11)*X11;
cog11_K_x = -X11*cos(phi11);
cog11_K_y = -X11*sin(phi11);

%Define additional centers of gravity for the calculation of the
%acceleration and rotation
cog45_D_x = (-X4*cos(phi3)*m4 + X5*cos(phi4)*m5)/(m4+m5);
cog45_D_y = (-X4*sin(phi3)*m4 + X5*sin(phi4)*m5)/(m4+m5);
cog89_H_x = (-X8*cos(phi8)*m8 + X9*cos(phi9)*m9)/(m8+m9);
cog89_H_y = (-X8*sin(phi8)*m8 + X9*sin(phi9)*m9)/(m8+m9);

%Additional distances
d_ED_x = -r5*cos(phi4);
d_ED_y = -r5*sin(phi4);
d_IH_x = -r9*cos(phi9);
d_IH_y = -r9*sin(phi9);

% 3D omega (dphi) and alpha (ddphi) vectors)
omega1 = [zeros(size(phi2)) zeros(size(phi2)) dphi1];
omega2 = [zeros(size(phi2)) zeros(size(phi2)) dphi2];
omega3 = [zeros(size(phi2)) zeros(size(phi2)) dphi3];
omega4 = [zeros(size(phi2)) zeros(size(phi2)) dphi4];
omega5 = [zeros(size(phi2)) zeros(size(phi2)) dphi5];
omega6 = [zeros(size(phi2)) zeros(size(phi2)) dphi6];
omega8 = [zeros(size(phi2)) zeros(size(phi2)) dphi8];
omega9 = [zeros(size(phi2)) zeros(size(phi2)) dphi9];
omega10 = [zeros(size(phi2)) zeros(size(phi2)) dphi10];
omega11 = [zeros(size(phi2)) zeros(size(phi2)) dphi11];
alpha1 = [zeros(size(phi2)) zeros(size(phi2)) ddphi1];
alpha2 = [zeros(size(phi2)) zeros(size(phi2)) ddphi2];
alpha3 = [zeros(size(phi2)) zeros(size(phi2)) ddphi3];
alpha4 = [zeros(size(phi2)) zeros(size(phi2)) ddphi4];
alpha5 = [zeros(size(phi2)) zeros(size(phi2)) ddphi5];
alpha6 = [zeros(size(phi2)) zeros(size(phi2)) ddphi6];
alpha8 = [zeros(size(phi2)) zeros(size(phi2)) ddphi8];
alpha9 = [zeros(size(phi2)) zeros(size(phi2)) ddphi9];
alpha10 = [zeros(size(phi2)) zeros(size(phi2)) ddphi10];
alpha11 = [zeros(size(phi2)) zeros(size(phi2)) ddphi11];

% 3D model vectors
cog2_A_vec = [cog2_A_x      cog2_A_y        zeros(size(phi2))];
AB_vec = [r2*cos(phi1)      r2*sin(phi1)    zeros(size(phi2))];
cog3_B_vec = [cog3_B_x      cog3_B_y        zeros(size(phi2))];
cog45_D_vec = [cog45_D_x    cog45_D_y       zeros(size(phi2))];
DE_vec = [r5*cos(phi4)      r5*sin(phi4)    zeros(size(phi2))];
cog6_E_vec = [cog6_E_x      cog6_E_y        zeros(size(phi2))];
cog7_G_vec = [cog7_G_x      cog7_G_y        zeros(size(phi2))];
cog89_H_vec = [cog89_H_x    cog89_H_y       zeros(size(phi2))];
HI_vec = [r9*cos(phi9)      r9*sin(phi9)    zeros(size(phi2))];
cog10_I_vec = [cog10_I_x    cog10_I_y       zeros(size(phi2))];
cog11_K_vec = [cog11_K_x    cog11_K_y       zeros(size(phi2))];

% acceleration vectors
acc_2 = cross(omega1,cross(omega1,cog2_A_vec)) + cross(alpha1,cog2_A_vec);
acc_B = cross(omega1,cross(omega1, AB_vec)) + cross(alpha1, AB_vec);
acc_3 = acc_B + cross(omega2, cross(omega2,  cog3_B_vec)) + cross(alpha2, cog3_B_vec);
acc_45 = cross(omega3, cross(omega3, cog45_D_vec)) + cross(alpha3, cog45_D_vec);
acc_E = cross(omega4, cross(omega4, DE_vec)) + cross(alpha4, DE_vec);
acc_6 = acc_E + cross(omega5, cross(omega5, cog6_E_vec)) + cross(alpha5, cog6_E_vec);
acc_7 = cross(omega6, cross(omega6, cog7_G_vec)) + cross(alpha6, cog7_G_vec);
acc_89 = cross(omega8, cross(omega8, cog89_H_vec)) + cross(alpha8, cog89_H_vec);
acc_I = cross(omega9, cross(omega9, HI_vec)) + cross(alpha9, HI_vec);
acc_10 = acc_I + cross(omega10, cross(omega10, cog10_I_vec)) + cross(alpha10, cog10_I_vec);
acc_11 = cross(omega11, cross(omega11, cog11_K_vec)) + cross(alpha11, cog11_K_vec);

acc_2x = acc_2(:,1);
acc_2y = acc_2(:,2);
acc_3x = acc_3(:,1);
acc_3y = acc_3(:,2);
acc_45x = acc_45(:,1);
acc_45y = acc_45(:,2);
acc_6x = acc_6(:,1);
acc_6y = acc_6(:,2);
acc_7x = acc_7(:,1);
acc_7y = acc_7(:,2);
acc_89x = acc_89(:,1);
acc_89y = acc_89(:,2);
acc_10x = acc_10(:,1);
acc_10y = acc_10(:,2);
acc_11x = acc_11(:,1);
acc_11y = acc_11(:,2);

% **********************
% *** force analysis ***
% **********************
% allocate matrices for force (F) and moment (M)
F_A_x = zeros(size(phi2));
F_A_y = zeros(size(phi2));
F_23_x = zeros(size(phi2));
F_23_y = zeros(size(phi2));
F_C_x = zeros(size(phi2));
F_34 = zeros(size(phi2));
F_38 = zeros(size(phi2));
F_D_x = zeros(size(phi2));
F_D_y = zeros(size(phi2));
F_56_x = zeros(size(phi2));
F_56_y = zeros(size(phi2));
F_67_x = zeros(size(phi2));
F_67_y = zeros(size(phi2));
F_G_x = zeros(size(phi2));
F_G_y = zeros(size(phi2));
F_H_x = zeros(size(phi2));
F_H_y = zeros(size(phi2));
F_910_x = zeros(size(phi2));
F_910_y = zeros(size(phi2));
F_1011_x = zeros(size(phi2));
F_1011_y = zeros(size(phi2));
F_K_x = zeros(size(phi2));
F_K_y = zeros(size(phi2));
M_A = zeros(size(phi2));

% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
  A = [ 1   0   -1  0   0   0           0           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   1   0  -1   0   0           0           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   1   0   1 sin(phi3(k)) -sin(phi8(k)) 0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   1   0 -cos(phi3(k)) cos(phi8(k)) 0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0 -sin(phi3(k)) 0           1   0  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0  cos(phi3(k)) 0           0   1   0  -1   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0           0           0   0   1   0  -1   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0           0           0   0   0   1   0  -1   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0           0           0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0           0           0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0          sin(phi8(k)) 0   0   0   0   0   0   0   0   1   0  -1   0   0   0   0   0   0;
        0   0   0   0   0   0         -cos(phi8(k)) 0   0   0   0   0   0   0   0   0   1   0  -1   0   0   0   0   0;
        0   0   0   0   0   0           0           0   0   0   0   0   0   0   0   0   0   1   0  -1   0   0   0   0;
        0   0   0   0   0   0           0           0   0   0   0   0   0   0   0   0   0   0   1   0  -1   0   0   0;
        0   0   0   0   0   0           0           0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0;
        0   0   0   0   0   0           0           0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0;
        cog2_A_y(k) -cog2_A_x(k) -cog2_B_y(k) cog2_B_x(k)   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1;
        0   0   cog3_B_y(k) -cog3_B_x(k) cog3_C_y(k) cog3_C_y(k)*sin(phi3(k))+cog3_C_x(k)*cos(phi3(k)) -cog3_C_y(k)*sin(phi8(k))-cog3_C_x(k)*cos(phi8(k))   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0  -r4(k)   0   0   0  -d_ED_y(k)   d_ED_x(k)   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0   0   0   0   cog6_E_y(k) -cog6_E_x(k) -cog6_F_y(k) cog6_F_x(k)   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0   0   0   0   0   0   cog7_F_y(k) -cog7_F_x(k) cog7_G_y(k) -cog7_G_x(k)   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0  -r8(k)   0   0   0   0   0   0   0   0   0   0  -d_IH_y(k)  d_IH_x(k)    0   0   0   0   0;
        0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   cog10_I_y(k) -cog10_I_x(k) -cog10_J_y(k) cog10_J_x(k)   0   0   0;
        0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   cog11_J_y(k) -cog11_J_x(k) cog11_K_y(k) -cog11_K_x(k)   0];
        
  B = [ m2*acc_2x(k);
        m2*acc_2y(k);
        m3*acc_3x(k);
        m3*acc_3y(k);
        (m4+m5)*acc_45x(k);
        m45*acc_45y(k);
        m6wing*acc_6x(k)-F_wing*cos(phi5(k));
        m6wing*acc_6y(k)-F_wing*sin(phi5(k));
        m7*acc_7x(k);
        m7*acc_7y(k);     
        m89*acc_89x(k);
        m89*acc_89y(k);
        m10wing*acc_10x(k)-F_wing*cos(phi10(k));
        m10wing*acc_10y(k)-F_wing*sin(phi10(k));
        m11*acc_11x(k);
        m11*acc_11y(k);
        J2*ddphi1(k);
        J3*ddphi2(k);
        J45*ddphi3(k);
        J6*ddphi5(k) - (Y6+Ywing)*F_wing;
        J7*ddphi6(k);
        J89*ddphi8(k);
        J10*ddphi10(k) - (Y10-Ywing)*F_wing;
        J11*ddphi11(k) ];
    
    
    
    x = A\B;
    
    % save results
    F_A_x(k) = x(1);
    F_A_y(k) = x(2);
    F_23_x(k) = x(3);
    F_23_y(k) = x(4);
    F_C_x(k) = x(5);
    F_34(k) = x(6);
    F_38(k) = x(7);
    F_D_x(k) = x(8);
    F_D_y(k) = x(9);
    F_56_x(k) = x(10);
    F_56_y(k) = x(11);
    F_67_x(k) = x(12);
    F_67_y(k) = x(13);
    F_G_x(k) = x(14);
    F_G_y(k) = x(15);
    F_H_x(k) = x(16);
    F_H_y(k) = x(17);
    F_910_x(k) = x(18);
    F_910_y(k) = x(19);
    F_1011_x(k) = x(20);
    F_1011_y(k) = x(21);
    F_K_x(k) = x(22);
    F_K_y(k) = x(23);
    M_A(k) = x(24);
    
end



% **********************
% *** plot figures ***
% **********************

if fig_dyn_4bar
    
    figure
    subplot(221)
    plot(F_A_x,F_A_y),grid
    xlabel('F_A_x [N]')
    ylabel('F_A_y [N]')
    axis tight
    subplot(222)
    plot(F_23_x,F_23_y),grid
    xlabel('F_23_x [N]')
    ylabel('F_Q23_y [N]')
    axis tight
    subplot(223)
    plot(F_34,F_38),grid
    xlabel('F_34 [N]')
    ylabel('F_38[N]')
    axis tight
    subplot(224)
    plot(F_C_x,zeros(size(phi2))),grid
    xlabel('F_C_x [N]')
    ylabel('F_C_y [N]')
    axis tight
    
        figure
    subplot(221)
    plot(F_D_x,F_D_y),grid
    xlabel('F_D_x [N]')
    ylabel('F_D_y [N]')
    axis tight
    subplot(222)
    plot(F_56_x,F_56_y),grid
    xlabel('F_56_x [N]')
    ylabel('F_56_y [N]')
    axis tight
    subplot(223)
    plot(F_67_x,F_67_y),grid
    xlabel('F_67_x [N]')
    ylabel('F_67_y [N]')
    axis tight
    subplot(224)
    plot(F_G_x,F_G_y),grid
    xlabel('F_G_x [N]')
    ylabel('F_G_y [N]')
    axis tight
    
    figure
    subplot(221)
    plot(F_H_x,F_H_y),grid
    xlabel('F_H_x [N]')
    ylabel('F_H_y [N]')
    axis tight
    subplot(222)
    plot(F_910_x,F_910_y),grid
    xlabel('F_910_x [N]')
    ylabel('F_910_y [N]')
    axis tight
    subplot(223)
    plot(F_1011_x,F_1011_y),grid
    xlabel('F_1011_x [N]')
    ylabel('F_1011_y [N]')
    axis tight
    subplot(224)
    plot(F_K_x,F_K_y),grid
    xlabel('F_K_x [N]')
    ylabel('F_K_y [N]')
    axis tight    
    
    figure
    plot(t,M_A)
    ylabel('M_A [N-m]')
    xlabel('t [s]')
    
end


