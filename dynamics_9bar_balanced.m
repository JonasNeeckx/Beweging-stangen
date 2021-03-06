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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F_A_x, F_A_y, F_23_x, F_23_y, F_C_x, F_34, F_38, F_D_x, F_D_y, F_56_x, F_56_y, F_67_x, F_67_y,...
    F_G_x, F_G_y, F_H_x, F_H_y, F_910_x, F_910_y, F_1011_x, F_1011_y, F_K_x, F_K_y, M_A] = ...
dynamics_9bar_balanced(phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,...
dphi2,dphi3,dphi4, dphi6, dphi7, dphi8,  dphi10, dphi11,...
ddphi2,ddphi3,ddphi4,ddphi6, ddphi7, ddphi8,  ddphi10, ddphi11,...
r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,rwing,rmax4, rmax8, ...
m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,mwing, ...
X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,Xwing, ...
Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Ywing, ...
J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,t,fig_dyn_9bar)


% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.

%As angle phi4 and phi9 are equal to respectivily phi3 and phi8 apart from
%a constant value, their dphi and ddphi is not given.
dphi5 = dphi4;
ddphi5 = ddphi4;
dphi9 = dphi8;
ddphi9 = ddphi8;

%Thrust on the wing (simplified to a force in the middle of the wing)
F_wing = 0;% 0.12;

% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P
cog2_A_x = zeros(size(phi2));
cog2_A_y = zeros(size(phi2));
cog2_B_x = -X2*cos(phi2);
cog2_B_y = -X2*sin(phi2);
cog3_B_x = -X3*cos(phi3);
cog3_B_y = -X3*sin(phi3);
cog3_C_x = X3*cos(phi3);
cog3_C_y = X3*sin(phi3); 
cog6_E_x = cos(phi6)*X6 - sin(phi6)*Y6;
cog6_E_y = sin(phi6)*X6 + cos(phi6)*Y6;
cog6_F_x = -cos(phi6)*X6 - sin(phi6)*Y6;
cog6_F_y = -sin(phi6)*X6 + cos(phi6)*Y6;
cog7_F_x = cos(phi7)*X7;
cog7_F_y = sin(phi7)*X7;
cog7_G_x = -cos(phi7)*X7;
cog7_G_y = -sin(phi7)*X7;
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
cog45_D_x = zeros(size(phi2));
cog45_D_y = zeros(size(phi2));
cog89_H_x = zeros(size(phi2));
cog89_H_y = zeros(size(phi2));

%Additional distances d_AB_i from point A to B over axis i
d_ED_x = -r5*cos(phi5);
d_ED_y = -r5*sin(phi5);
d_IH_x = -r9*cos(phi9);
d_IH_y = -r9*sin(phi9);

%The masses and moments of inertia of bonded links are added
%together
m45 = m4 + m5;
m89 = m8 + m9;
m6wing = m6 + mwing;
m10wing = m10 + mwing;
J45 = J4 + J5;
J45cog = J45 + m45*(cog45_D_x.^2+cog45_D_y.^2);
J45cog = J45cog(1);
J89 = J8 + J9;
J89cog = J89 + m89*(cog89_H_x.^2+cog89_H_y.^2);
J89cog = J89cog(1);

% 3D omega (dphi) and alpha (ddphi) vectors)
omega2 = [zeros(size(phi2)) zeros(size(phi2)) dphi2];
omega3 = [zeros(size(phi2)) zeros(size(phi2)) dphi3];
omega4 = [zeros(size(phi2)) zeros(size(phi2)) dphi4];
omega5 = [zeros(size(phi2)) zeros(size(phi2)) dphi5];
omega6 = [zeros(size(phi2)) zeros(size(phi2)) dphi6];
omega7 = [zeros(size(phi2)) zeros(size(phi2)) dphi7];
omega8 = [zeros(size(phi2)) zeros(size(phi2)) dphi8];
omega9 = [zeros(size(phi2)) zeros(size(phi2)) dphi9];
omega10 = [zeros(size(phi2)) zeros(size(phi2)) dphi10];
omega11 = [zeros(size(phi2)) zeros(size(phi2)) dphi11];
alpha2 = [zeros(size(phi2)) zeros(size(phi2)) ddphi2];
alpha3 = [zeros(size(phi2)) zeros(size(phi2)) ddphi3];
alpha4 = [zeros(size(phi2)) zeros(size(phi2)) ddphi4];
alpha5 = [zeros(size(phi2)) zeros(size(phi2)) ddphi5];
alpha6 = [zeros(size(phi2)) zeros(size(phi2)) ddphi6];
alpha7 = [zeros(size(phi2)) zeros(size(phi2)) ddphi7];
alpha8 = [zeros(size(phi2)) zeros(size(phi2)) ddphi8];
alpha9 = [zeros(size(phi2)) zeros(size(phi2)) ddphi9];
alpha10 = [zeros(size(phi2)) zeros(size(phi2)) ddphi10];
alpha11 = [zeros(size(phi2)) zeros(size(phi2)) ddphi11];

% 3D model vectors
cog2_A_vec = [cog2_A_x      cog2_A_y        zeros(size(phi2))];
AB_vec = [r2*cos(phi2)      r2*sin(phi2)    zeros(size(phi2))];
cog3_B_vec = [cog3_B_x      cog3_B_y        zeros(size(phi2))];
cog45_D_vec = [cog45_D_x    cog45_D_y       zeros(size(phi2))];
DE_vec = [r5*cos(phi5)      r5*sin(phi5)    zeros(size(phi2))];
cog6_E_vec = [cog6_E_x      cog6_E_y        zeros(size(phi2))];
cog7_G_vec = [cog7_G_x      cog7_G_y        zeros(size(phi2))];
cog89_H_vec = [cog89_H_x    cog89_H_y       zeros(size(phi2))];
HI_vec = [r9*cos(phi9)      r9*sin(phi9)    zeros(size(phi2))];
cog10_I_vec = [cog10_I_x    cog10_I_y       zeros(size(phi2))];
cog11_K_vec = [cog11_K_x    cog11_K_y       zeros(size(phi2))];

% acceleration vectors
acc_2 = cross(omega2,cross(omega2,cog2_A_vec)) + cross(alpha2,cog2_A_vec);
acc_B = cross(omega2,cross(omega2, AB_vec)) + cross(alpha2, AB_vec);
acc_3 = acc_B + cross(omega3, cross(omega3,  cog3_B_vec)) + cross(alpha3, cog3_B_vec);
acc_45 = cross(omega4, cross(omega4, cog45_D_vec)) + cross(alpha4, cog45_D_vec);
acc_E = cross(omega5, cross(omega5, DE_vec)) + cross(alpha5, DE_vec);
acc_6 = acc_E + cross(omega6, cross(omega6, cog6_E_vec)) + cross(alpha6, cog6_E_vec);
acc_7 = cross(omega7, cross(omega7, cog7_G_vec)) + cross(alpha7, cog7_G_vec);
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
        0   0   1   0   1 sin(phi4(k)) -sin(phi8(k)) 0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   1   0 -cos(phi4(k)) cos(phi8(k)) 0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0 -sin(phi4(k)) 0           1   0  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0  cos(phi4(k)) 0           0   1   0  -1   0   0   0   0   0   0   0   0   0   0   0   0   0;
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
        0   0   cog3_B_y(k) -cog3_B_x(k) cog3_C_y(k) cog3_C_y(k)*sin(phi4(k))+cog3_C_x(k)*cos(phi4(k)) -cog3_C_y(k)*sin(phi8(k))-cog3_C_x(k)*cos(phi8(k))  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0  -r4(k)   0   0   0  -d_ED_y(k)   d_ED_x(k)   0   0   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0   0   0   0   cog6_E_y(k) -cog6_E_x(k) -cog6_F_y(k) cog6_F_x(k)   0   0   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0   0   0   0   0   0   cog7_F_y(k) -cog7_F_x(k) cog7_G_y(k) -cog7_G_x(k)   0   0   0   0   0   0   0   0   0;
        0   0   0   0   0   0   r8(k)   0   0   0   0   0   0   0   0   0   0  -d_IH_y(k)  d_IH_x(k)    0   0   0   0   0;
        0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   cog10_I_y(k) -cog10_I_x(k) -cog10_J_y(k) cog10_J_x(k)   0   0   0;
        0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   cog11_J_y(k) -cog11_J_x(k) cog11_K_y(k) -cog11_K_x(k)   0];
        
  B = [ m2*acc_2x(k);
        m2*acc_2y(k);
        m3*acc_3x(k);
        m3*acc_3y(k);
        m45*acc_45x(k);
        m45*acc_45y(k);
        m6wing*acc_6x(k)-F_wing*cos(phi6(k));
        m6wing*acc_6y(k)-F_wing*sin(phi6(k));
        m7*acc_7x(k);
        m7*acc_7y(k);     
        m89*acc_89x(k);
        m89*acc_89y(k);
        m10wing*acc_10x(k)-F_wing*cos(phi10(k));
        m10wing*acc_10y(k)-F_wing*sin(phi10(k));
        m11*acc_11x(k);
        m11*acc_11y(k);
        J2*ddphi2(k);
        J3*ddphi3(k);
        J45*ddphi4(k);
        J6*ddphi6(k) - (Y6+Ywing)*F_wing;
        J7*ddphi7(k);
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

%**********************
%***** controle *******
%**********************

% velocity vectors
vel_2 = cross(omega2,cog2_A_vec);
vel_B = cross(omega2, AB_vec);
vel_3 = vel_B + cross(omega3,  cog3_B_vec);
vel_45 = cross(omega4, cog45_D_vec);
vel_E = cross(omega5, DE_vec);
vel_6 = vel_E + cross(omega6, cog6_E_vec);
vel_7 = cross(omega7, cog7_G_vec);
vel_89 = cross(omega8, cog89_H_vec);
vel_I = cross(omega9, HI_vec);
vel_10 = vel_I + cross(omega10, cog10_I_vec);
vel_11 = cross(omega11, cog11_K_vec);

vel_2x = vel_2(:,1);
vel_2y = vel_2(:,2);
vel_3x = vel_3(:,1);
vel_3y = vel_3(:,2);
vel_45x = vel_45(:,1);
vel_45y = vel_45(:,2);
vel_6x = vel_6(:,1);
vel_6y = vel_6(:,2);
vel_7x = vel_7(:,1);
vel_7y = vel_7(:,2);
vel_89x = vel_89(:,1);
vel_89y = vel_89(:,2);
vel_10x = vel_10(:,1);
vel_10y = vel_10(:,2);
vel_11x = vel_11(:,1);
vel_11y = vel_11(:,2);

M_A_Energy = zeros(size(phi3));
F_x_Work = zeros(size(phi4));
F_y_Work = zeros(size(phi3));

for k = 1:t_size
M_A_Energy(k) = (m2 * (acc_2x(k)*vel_2x(k) + acc_2y(k)*vel_2y(k))+m3 * (acc_3x(k)*vel_3x(k) + acc_3y(k)*vel_3y(k))+...
    m45 * (acc_45x(k)*vel_45x(k) + acc_45y(k)*vel_45y(k))+m6wing * (acc_6x(k)*vel_6x(k) + acc_6y(k)*vel_6y(k))+...
    m7 * (acc_7x(k)*vel_7x(k) + acc_7y(k)*vel_7y(k))+m89 * (acc_89x(k)*vel_89x(k) + acc_89y(k)*vel_89y(k))+...
    m10wing * (acc_10x(k)*vel_10x(k) + acc_10y(k)*vel_10y(k))+m11 * (acc_11x(k)*vel_11x(k) + acc_11y(k)*vel_11y(k))+...
    J2*dphi2(k)*ddphi2(k)+J3*dphi3(k)*ddphi3(k)+J45cog*dphi4(k)*ddphi4(k)+J6*dphi6(k)*ddphi6(k)+J7*dphi7(k)*ddphi7(k)+...
    J89cog*dphi8(k)*ddphi8(k)+J10*dphi10(k)*ddphi10(k)+J11*dphi11(k)*ddphi11(k))/dphi2(k); 

F_x_Work(k) = m2*acc_2x(k) + m3*acc_3x(k) + m45*acc_45x(k) + m6wing*acc_6x(k) + ...
    m7*acc_7x(k) + m89*acc_89x(k) + m10wing*acc_10x(k) + m11*acc_11x(k);

F_y_Work(k) = m2*acc_2y(k) + m3*acc_3y(k) + m45*acc_45y(k) + m6wing*acc_6y(k) + ...
    m7*acc_7y(k) + m89*acc_89y(k) + m10wing*acc_10y(k) + m11*acc_11y(k);
end

% **********************
% *** plot figures ***
% **********************

if fig_dyn_9bar
    
    figure
    subplot(211)
    plot(F_A_x,F_A_y),grid
    xlabel('F_A_x_b_a_l_a_n_c_e_d [N]')
    ylabel('F_A_y_b_a_l_a_n_c_e_d [N]')
    axis tight
    
    
    subplot(212)
    plot(t,M_A)
    ylabel('M_A_b_a_l_a_n_c_e_d [N-m]')
    xlabel('t [s]')  
    
    
end


