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

function [phi2,phi3,phi5,phi6,phi8,phi10,phi11,dphi2,dphi3,dphi5,dphi6,dphi8,dphi10,dphi11,ddphi2,ddphi3,ddphi5,ddphi6,ddphi8,ddphi10,ddphi11,r8,dr8,ddr8,r13,dr13,ddr13,r4,dr4,ddr4] = kinematics_4bar(r2,r3,r5,r6,r7,r9,r10,r11,r14x,r14y,r47y,r18x,r18y,r811y,phi1,dphi1,ddphi1,phi2_init,phi3_init,phi5_init,phi6_init,phi8_init,phi10_init,phi11_init,r13_init,r4_init,r8_init,t,fig_kin_4bar);

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
phi2 = zeros(size(t));
phi3 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
phi8 = zeros(size(t));
phi10 = zeros(size(t));
phi11 = zeros(size(t));
dphi2 = zeros(size(t));
dphi3 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dphi8 = zeros(size(t));
dphi10 = zeros(size(t));
dphi11 = zeros(size(t));
ddphi2 = zeros(size(t));
ddphi3 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi8 = zeros(size(t));
ddphi10 = zeros(size(t));
ddphi11 = zeros(size(t));
r13 = zeros(size(t));
dr13 = zeros(size(t));
ddr13 = zeros(size(t));
r4 = zeros(size(t));
dr4 = zeros(size(t));
ddr4 = zeros(size(t));
r8 = zeros(size(t));
dr8 = zeros(size(t));
ddr8 = zeros(size(t));

phi4 = zeros(size(t));
phi9 = zeros(size(t));
% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off','TolFun',1e-11,'TolX',eps,'MaxIter',400);

% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles phi3 and phi4
    % argument optim options: parameters for fsolve
    % argument phi2(k): input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
    % argument a1 ... phi1: constants
    % return value x: solution for the unknown angles phi3 and phi4
    % return exitflag: indicates convergence of algorithm
    gamma = 330*pi/180;
    [x, fval, exitflag]=fsolve('loop_closure_eqs',[phi2_init phi3_init phi5_init phi6_init phi8_init phi10_init phi11_init r13_init r4_init r8_init]',optim_options,phi1(k),r2,r3,r5,r6,r7,r9,r10,r11,gamma,r14x,r14y,r47y,r18x,r18y,r811y);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % save results of fsolve
    phi2(k)=x(1);
    phi3(k)=x(2);
    phi5(k)=x(3);
    phi6(k)=x(4);
    phi8(k)=x(5);
    phi10(k)=x(6);
    phi11(k)=x(7);
    r13(k)=x(8);
    r4(k)=x(9);
    r8(k)=x(10);
    
    phi4(k) = phi3(k) + gamma;
    phi9(k) = phi8(k) + 2*pi - gamma;
    % *** velocity analysis ***
    
    A = [-r3*sin(phi2(k)),0,0,0,0,0,0,0,0,0;
        r3*cos(phi2(k)),0,0,0,0,0,0,1,0,0;
        0,r4(k)*sin(phi3(k)),0,0,0,0,0,0,-cos(phi3(k)),0;
        0,-r4(k)*cos(phi3(k)),0,0,0,0,0,1,-sin(phi3(k)),0;
        0,0,0,0,r8(k)*sin(phi8(k)),0,0,0,0,-cos(phi8(k));
        0,0,0,0,-r8(k)*cos(phi8(k)),0,0,1,0,-sin(phi8(k));
        0,-r5*sin(phi4(k)),-r6*sin(phi5(k)),-r7*sin(phi6(k)),0,0,0,0,0,0;
        0,r5*cos(phi4(k)),r6*cos(phi5(k)),r7*cos(phi6(k)),0,0,0,0,0,0;
        0,0,0,0,-r9*sin(phi9(k)),-r10*sin(phi10(k)),-r11*sin(phi11(k)),0,0,0;
        0,0,0,0,r9*cos(phi9(k)),r10*cos(phi10(k)),r11*cos(phi11(k)),0,0,0];
        
    B = [r2*sin(phi1(k))*dphi1(k);
        -r2*cos(phi1(k))*dphi1(k);
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0];
     
    x = A\B;
    % save results
    dphi2(k)=x(1);
    dphi3(k)=x(2);
    dphi5(k)=x(3);
    dphi6(k)=x(4);
    dphi8(k)=x(5);
    dphi10(k)=x(6);
    dphi11(k)=x(7);
    dr13(k)=x(8);
    dr4(k)=x(9);
    dr8(k)=x(10);
    
    
    % *** acceleration analysis ***
    
        
      A = A;
    
      C = [r2*cos(phi1(k))*dphi1(k)^2+r2*sin(phi1(k))*ddphi1(k)+r3*cos(phi2(k))*dphi2(k)^2;
          -r2*cos(phi1(k))*ddphi1(k)+r2*sin(phi1(k))*dphi1(k)^2+r3*sin(phi2(k))*dphi2(k)^2;
          -dphi3(k)^2*r4(k)*cos(phi3(k))-dphi3(k)*sin(phi3(k))*dr4(k)-dr4(k)*dphi3(k)*sin(phi3(k));
          dphi3(k)^2*r4(k)*sin(phi3(k))+dphi3(k)*cos(phi3(k))*dr4(k)+dr4(k)*dphi3(k)*cos(phi3(k));
          -dphi8(k)^2*r8(k)*cos(phi8(k))-dphi8(k)*sin(phi8(k))*dr8(k)-dr8(k)*dphi8(k)*sin(phi8(k));
          dphi8(k)^2*r8(k)*sin(phi8(k))+dphi8(k)*cos(phi8(k))*dr8(k)+dr8(k)*dphi8(k)*cos(phi8(k));
          r5*cos(phi4(k))*dphi3(k)^2+r6*cos(phi5(k))*dphi5(k)^2+r7*cos(phi6(k))*dphi6(k)^2;
          r5*sin(phi4(k))*dphi3(k)^2+r6*sin(phi5(k))*dphi5(k)^2+r7*sin(phi6(k))*dphi6(k)^2;
          r9*cos(phi9(k))*dphi8(k)^2+r10*cos(phi10(k))*dphi10(k)^2+r11*cos(phi11(k))*dphi11(k)^2;
          r9*sin(phi9(k))*dphi8(k)^2+r10*sin(phi10(k))*dphi10(k)^2+r11*cos(phi11(k))*dphi11(k)^2
          ];
    
    x = A\C;
    % save results
    ddphi2(k)=x(1);
    ddphi3(k)=x(2);
    ddphi5(k)=x(3);
    ddphi6(k)=x(4);
    ddphi8(k)=x(5);
    ddphi10(k)=x(6);
    ddphi11(k)=x(7);
    ddr13(k)=x(8);
    ddr4(k)=x(9);
    ddr8(k)=x(10);
    
    
    % *** calculate initial values for next iteration step ***
    phi2_init = phi2(k)+Ts*dphi2(k);
    phi3_init = phi3(k)+Ts*dphi3(k);
    phi5_init = phi5(k)+Ts*dphi5(k);
    phi6_init = phi6(k)+Ts*dphi6(k);
    phi8_init = phi8(k)+Ts*dphi8(k);
    phi10_init = phi10(k)+Ts*dphi10(k);
    phi11_init = phi11(k)+Ts*dphi11(k);
    r13_init = r13(k)+Ts*dr13(k);
    r4_init = r4(k)+Ts*dr4(k);
    r8_init = r8(k)+Ts*dr8(k);
    
    
end % loop over positions



% *** create movie ***

% point P = fixed
A = 0;
% point S = fixed
D = A -0.007746*exp(j*63.14*pi/180);
G = A -0.010983*exp(j*71.4166*pi/180);
H = A -0.007746*exp(j*63.14*pi/180);
K = A -0.010983*exp(j*71.4166*pi/180);


% define which positions we want as frames in our movie
frames = 40;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -0.025;
y_bottom = -0.025;
x_right = 0.025;
y_top = 0.025;

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    B = A + r2 * exp(j*phi1(index));
    C = B - r3 * exp(j*phi2(index));
    
    loop1 = [A B C A];
    
    D = C + r4(index) * exp(j*phi3(index));
    
    loop2 = [A C D A]
    
    H = C + r8(index) * exp(j*phi8(index));
    
    loop3 = [A C H A]
    
    E = D - r5*exp(j*phi4(index));
    F = E + r6*exp(j*phi5(index));
    G = F - r7*exp(j*phi6(index));
    
    loop4 = [D E F G D]
    
    I = H -r9*exp(j*phi9(index));
    J = I + r10*exp(j*phi10(index));
    K = J - r11*exp(j*phi11(index));
    
    loop5 = [H I J K H]
    figure(10)
    clf
    hold on
    plot(real(loop1),imag(loop1),'-o')
    plot(real(loop2),imag(loop2),'-o')
    plot(real(loop3),imag(loop3),'-o')
    plot(real(loop4),imag(loop4),'-o')
    plot(real(loop5),imag(loop5),'-o')
    
    axis(movie_axes);
    xlabel('x [m]')
    ylabel('y [m]')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on')% set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save fourbar_movie Movie
close(10)


% % *** plot figures ***
%     
%     figure
%     subplot(311)
%     plot(t,phi2)
%     ylabel('\phi_2 [rad]')
%     subplot(312)
%     plot(t,phi3)
%     ylabel('\phi_3 [rad]')
%     subplot(313)
%     plot(t,phi4)
%     ylabel('\phi_4 [rad]')
%     xlabel('t [s]')
%     
%     figure
%     subplot(311)
%     plot(t,dphi2)
%     ylabel('d\phi_2 [rad/s]')
%     subplot(312)
%     plot(t,dphi3)
%     ylabel('d\phi_3 [rad/s]')
%     subplot(313)
%     plot(t,dphi4)
%     ylabel('d\phi_4 [rad/s]')
%     xlabel('t [s]')
%     
%     figure
%     subplot(311)
%     plot(t,ddphi2)
%     ylabel('dd\phi_2 [rad/s^2]')
%     subplot(312)
%     plot(t,ddphi3)
%     ylabel('dd\phi_3 [rad/s^2]')
%     subplot(313)
%     plot(t,ddphi4)
%     ylabel('dd\phi_4 [rad/s^2]')
%     xlabel('t [s]')
end



