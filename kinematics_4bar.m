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

function [phi2,phi3,phi4,phi5,phi6,phi8,phi9,phi10,phi11,dphi2,dphi3,dphi4,dphi5,dphi6,dphi8,dphi9,dphi10,dphi11,ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi8,ddphi9,ddphi10,ddphi11] = kinematics_4bar(r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r14x,r14y,r18x,r18y,r811y,phi1,dphi1,ddphi1,phi2_init,phi3_init,phi4_init,phi5_init,phi6_init,phi8_init,phi9_init,phi10_init,phi11_init,r8_init,t,fig_kin_4bar)

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
phi2 = zeros(size(t));
phi3 = zeros(size(t));
phi4 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
phi8 = zeros(size(t));
phi9 = zeros(size(t));
phi10 = zeros(size(t));
phi11 = zeros(size(t));
dphi2 = zeros(size(t));
dphi3 = zeros(size(t));
dphi4 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dphi8 = zeros(size(t));
dphi9 = zeros(size(t));
dphi10 = zeros(size(t));
dphi11 = zeros(size(t));
ddphi2 = zeros(size(t));
ddphi3 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi8 = zeros(size(t));
ddphi9 = zeros(size(t));
ddphi10 = zeros(size(t));
ddphi11 = zeros(size(t));

% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off');

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
    [x, fval, exitflag]=fsolve('loop_closure_eqs',[phi2_init phi3_init phi4_init phi5_init phi6_init phi8_init phi9_init phi10_init phi11_init r8_init]',optim_options,phi1(k),r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r14x,r14y,r18x,r18y,r811y);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % save results of fsolve
    phi2(k)=x(1);
    phi3(k)=x(2);
    phi5(k)=x(4);
    phi6(k)=x(5);
    phi8(k)=x(6);
    phi10(k)=x(8);
    phi11(k)=x(9);
    
    
    
    % *** velocity analysis ***
    
    A = 
    B = 
     
    x = A\B;
    
    % save results
    dphi2(k)=x(1);
    dphi3(k)=x(2);
    dphi4(k)=x(3);
    dphi5(k)=x(4);
    dphi6(k)=x(5);
    dphi8(k)=x(6);
    dphi9(k)=x(7);
    dphi10(k)=x(8);
    dphi11(k)=x(9);
    
    
    % *** acceleration analysis ***
    
    A =
    B = 
    x = A\B;
    % save results
    ddphi2(k) = x(1);
    ddphi3(k) = x(2);
    ddphi4(k) = x(3);
    ddphi5(k) = x(4);
    ddphi6(k) = x(5);
    ddphi8(k) = x(6);
    ddphi9(k) = x(7);
    ddphi10(k) = x(8);
    ddphi11(k) = x(9);
    
    
    % *** calculate initial values for next iteration step ***
    phi2_init = phi2(k)+Ts*dphi2(k);
    phi3_init = phi3(k)+Ts*dphi3(k);
    phi4_init = phi4(k)+Ts*dphi4(k);
    phi5_init = phi5(k)+Ts*dphi5(k);
    phi6_init = phi6(k)+Ts*dphi6(k);
    phi8_init = phi8(k)+Ts*dphi8(k);
    phi9_init = phi9(k)+Ts*dphi9(k);
    phi10_init = phi10(k)+Ts*dphi10(k);
    phi11_init = phi11(k)+Ts*dphi11(k);
    
    
end % loop over positions



% *** create movie ***

% point P = fixed
P = 0;
% point S = fixed
S = r1*exp(j*phi1);
% define which positions we want as frames in our movie
frames = 40;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -1.5*r2;
y_bottom = -1.5*max(r2,r4);
x_right = r1+1.5*r4;
y_top = 1.5*max(r2,r4);

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    Q = P + r2 * exp(j*phi2(index));
    R1 = Q + r3 * exp(j*phi3(index));
    R2 = S + r4 * exp(j*phi4(index));
    
    loop1 = [P Q R1 R2 S];
    
    figure(10)
    clf
    hold on
    plot(real(loop1),imag(loop1),'-o')
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save fourbar_movie Movie
close(10)


% *** plot figures ***

if fig_kin_4bar
    
    %plot assembly at a certain timestep 
    index = 1; %select 1st timestep
    P = 0;
    S = r1*exp(j*phi1);
    Q = P + r2 * exp(j*phi2(index));
    R = Q + r3 * exp(j*phi3(index));
    
    figure
    assembly=[P, Q, R, S];
    plot(real(assembly),imag(assembly),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    figure
    subplot(311)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(312)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(313)
    plot(t,phi4)
    ylabel('\phi_4 [rad]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    subplot(312)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    subplot(313)
    plot(t,dphi4)
    ylabel('d\phi_4 [rad/s]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(312)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(313)
    plot(t,ddphi4)
    ylabel('dd\phi_4 [rad/s^2]')
    xlabel('t [s]')
end



