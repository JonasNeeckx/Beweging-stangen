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

function [phi2,phi3,phi4,phi5,phi6,phi8,phi9,phi10,phi11,dphi2,dphi3,dphi5,dphi6,dphi8,dphi10,dphi11,ddphi2,ddphi3,ddphi5,ddphi6,ddphi8,ddphi10,ddphi11,r8,dr8,ddr8,r13,dr13,ddr13,r4,dr4,ddr4] = kinematics_4bar(r2,r3,r5,r6,r7,r9,r10,r11,r14x,r14y,r47y,r18x,r18y,r811y,phi1,dphi1,ddphi1,phi2_init,phi3_init,phi5_init,phi6_init,phi8_init,phi10_init,phi11_init,r13_init,r4_init,r8_init,t,fig_kin_4bar);

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
    
    A = [r3*sin(phi2(k)),0,0,0,0,0,0,0,0,0;
        -r3*cos(phi2(k)),0,0,0,0,0,0,1,0,0;
        0,r4(k)*sin(phi3(k)),0,0,0,0,0,0,-cos(phi3(k)),0;
        0,-r4(k)*cos(phi3(k)),0,0,0,0,0,1,-sin(phi3(k)),0;
        0,0,0,0,r8(k)*sin(phi8(k)),0,0,0,0,-cos(phi8(k));
        0,0,0,0,-r8(k)*cos(phi8(k)),0,0,1,0,-sin(phi8(k));
        0,-r5*sin(phi4(k)),-r6*sin(phi5(k)),-r7*sin(phi6(k)),0,0,0,0,0,0;
        0,r5*cos(phi4(k)),r6*cos(phi5(k)),r7*cos(phi6(k)),0,0,0,0,0,0;
        0,0,0,0,-r9*sin(phi9(k)),-r10*sin(phi10(k)),-r11*sin(phi11(k)),0,0,0;
        0,0,0,0,r9*cos(phi9(k)),r10*cos(phi10(k)),r11*cos(phi11(k)),0,0,0];
        
    B = [+r2*sin(phi1(k))*dphi1(k);
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
    
        
      A = [r3*sin(phi2(k)),0,0,0,0,0,0,0,0,0;
          -r3*cos(phi2(k)),0,0,0,0,0,0,1,0,0;
          0,r4(k)*sin(phi3(k)),0,0,0,0,0,0,-cos(phi3(k)),0;
          0,-r4(k)*cos(phi3(k)),0,0,0,0,0,1,-sin(phi3(k)),0;
          0,0,0,0,r8(k)*sin(phi8(k)),0,0,0,0,-cos(phi8(k));
          0,0,0,0,-r8(k)*cos(phi8(k)),0,0,1,0,-sin(phi8(k));
          0,-r5*sin(phi4(k)),-r6*sin(phi5(k)),-r7*sin(phi6(k)),0,0,0,0,0,0;
          0,r5*cos(phi4(k)),r6*cos(phi5(k)),r7*cos(phi6(k)),0,0,0,0,0,0;
          0,0,0,0,-r9*sin(phi9(k)),-r10*sin(phi10(k)),-r11*sin(phi11(k)),0,0,0;
          0,0,0,0,r9*cos(phi9(k)),r10*cos(phi10(k)),r11*cos(phi11(k)),0,0,0];
    
      B = [r2*cos(phi1(k))*dphi1(k)^2+r2*sin(phi1(k))*ddphi1(k)-r3*cos(phi2(k))*dphi2(k)^2;
          -r2*cos(phi1(k))*ddphi1(k)+r2*sin(phi1(k))*dphi1(k)^2-r3*sin(phi2(k))*dphi2(k)^2;
          -dphi3(k)^2*r4(k)*cos(phi3(k))-dphi3(k)*sin(phi3(k))*dr4(k)-dr4(k)*dphi3(k)*sin(phi3(k));
          -dphi3(k)^2*r4(k)*sin(phi3(k))+dphi3(k)*cos(phi3(k))*dr4(k)+dr4(k)*dphi3(k)*cos(phi3(k));
          -dphi8(k)^2*r8(k)*cos(phi8(k))-dphi8(k)*sin(phi8(k))*dr8(k)-dr8(k)*dphi8(k)*sin(phi8(k));
          -dphi8(k)^2*r8(k)*sin(phi8(k))+dphi8(k)*cos(phi8(k))*dr8(k)+dr8(k)*dphi8(k)*cos(phi8(k));
          r5*cos(phi4(k))*dphi3(k)^2+r6*cos(phi5(k))*dphi5(k)^2+r7*cos(phi6(k))*dphi6(k)^2;
          r5*sin(phi4(k))*dphi3(k)^2+r6*sin(phi5(k))*dphi5(k)^2+r7*sin(phi6(k))*dphi6(k)^2;
          r9*cos(phi9(k))*dphi8(k)^2+r10*cos(phi10(k))*dphi10(k)^2+r11*cos(phi11(k))*dphi11(k)^2;
          r9*sin(phi9(k))*dphi8(k)^2+r10*sin(phi10(k))*dphi10(k)^2+r11*sin(phi11(k))*dphi11(k)^2
          ];
    
    x = A\B;
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

<<<<<<< HEAD


% *** create movie ***

% point P = fixed
A = 0;
% point S = fixed
D = A +sqrt(0.0035^2+0.00691^2)*exp(j*+(acosd(0.00691/(sqrt(0.0035^2+0.00691^2)))-90)*pi/180);
G = A +sqrt(0.0035^2+(0.00691+0.0035)^2)*exp(j*+(acosd((0.00691+0.0035)/sqrt(0.0035^2+(0.00691+0.0035)^2))-90)*pi/180);
H = A -sqrt(0.0035^2+0.00691^2)*exp(j*-(acosd(0.00691/sqrt(0.0035^2+0.00691^2))-90)*pi/180);
K = A -sqrt(0.0035^2+(0.00691+0.0035)^2)*exp(j*-(acosd((0.00691+0.0035)/sqrt(0.0035^2+(0.00691+0.0035)^2))-90)*pi/180);


% define which positions we want as frames in our movie
frames = 200;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -0.015;
y_bottom = -0.015;
x_right = 0.015;
y_top = 0.005;

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
    
    E = D + r5*exp(j*phi4(index));
    F = G + r7*exp(j*(phi6(index)-pi));

    I = H +r9*exp(j*phi9(index));
    J = K + r11*exp(j*(phi11(index)-pi));
    loop = [A B C D E F G F E D C H I J K];
    figure(10)
    clf
    hold on
    plot(real(loop),imag(loop),'-o')
    
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
if fig_kin_4bar
    
    %plot assembly at a certain timestep 
    B = A + r2 * exp(j*phi1(1));
    C = B - r3 * exp(j*phi2(1));
    
    E = D + r5*exp(j*phi4(1));
    F = G + r7*exp(j*(phi6(1)-pi));

    I = H +r9*exp(j*phi9(1));
    J = K + r11*exp(j*(phi11(1)-pi));%select 1st timestep
    
    figure
    assembly=[A B C D E F G F E D C H I J K];
    plot(real(assembly),imag(assembly),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
end
% *** plot figures ***
    figure('Name','Position: \phi_2 - \r13')
    subplot(4,3,1)
    plot(t,phi1)
    ylabel('\phi_1 [rad]')
    subplot(4,3,2)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(4,3,3)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(4,3,4)
    plot(t,phi5)
    ylabel('\phi_5 [rad]')
    xlabel('t [s]')
    subplot(4,3,5)
    plot(t,phi6)
    ylabel('\phi_6 [rad]')
    xlabel('t [s]')
    subplot(4,3,6)
    plot(t,phi8)
    ylabel('\phi_8 [rad]')
    xlabel('t [s]')
    subplot(4,3,7)
    plot(t,phi10)
    ylabel('\phi_10 [rad]')
    xlabel('t [s]')
    subplot(4,3,8)
    plot(t,phi11)
    ylabel('\phi_11 [rad]')
    xlabel('t [s]')
    subplot(4,3,9)
    plot(t,r4)
    ylabel('r4 [m]')
    xlabel('t [s]')
    subplot(4,3,10)
    plot(t,r8)
    ylabel('r8 [m]')
    xlabel('t [s]')
    subplot(4,3,11)
    plot(t,r13)
    ylabel('r13 [m]')
    xlabel('t [s]')
    
    figure
    subplot(4,3,1)
    plot(t,dphi1)
    ylabel('d\phi_1 [rad/s]')
    xlabel('t [s]')
    subplot(4,3,2)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    xlabel('t [s]')
    subplot(4,3,3)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    xlabel('t [s]')
    subplot(4,3,4)
    plot(t,dphi5)
    ylabel('d\phi_5 [rad/s]')
    xlabel('t [s]')
    subplot(4,3,5)
    plot(t,dphi6)
    ylabel('d\phi_6 [rad/s]')
    xlabel('t [s]')
    subplot(4,3,6)
    plot(t,dphi8)
    ylabel('d\phi_8 [rad/s]')
    xlabel('t [s]')
    subplot(4,3,7)
    plot(t,dphi10)
    ylabel('d\phi_10 [rad/s]')
    xlabel('t [s]')
    subplot(4,3,8)
    plot(t,dphi11)
    ylabel('d\phi_11 [rad/s]')
    xlabel('t [s]')
    subplot(4,3,9)
    plot(t,dr4)
    ylabel('dr4 [m/s]')
    xlabel('t [s]')
    subplot(4,3,10)
    plot(t,dr8)
    ylabel('dr8 [m/s]')
    xlabel('t [s]')
    subplot(4,3,11)
    plot(t,dr13)
    ylabel('dr13 [m/s]')
    xlabel('t [s]')
    
    figure
    subplot(4,3,1)
    plot(t,ddphi1)
    ylabel('dd\phi_1 [rad/s^2]')
    xlabel('t [s]')
    subplot(4,3,2)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    xlabel('t [s]')
    subplot(4,3,3)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    xlabel('t [s]')
    subplot(4,3,4)
    plot(t,ddphi5)
    ylabel('dd\phi_5 [rad/s^2]')
    xlabel('t [s]')
    subplot(4,3,5)
    plot(t,ddphi6)
    ylabel('dd\phi_6 [rad/s^2]')
    xlabel('t [s]')
    subplot(4,3,6)
    plot(t,ddphi8)
    ylabel('dd\phi_8 [rad/s^2]')
    xlabel('t [s]')
    subplot(4,3,7)
    plot(t,ddphi10)
    ylabel('dd\phi_10 [rad/s^2]')
    xlabel('t [s]')
    subplot(4,3,8)
    plot(t,ddphi11)
    ylabel('dd\phi_11 [rad/s^2]')
    xlabel('t [s]')
    subplot(4,3,9)
    plot(t,ddr4)
    ylabel('ddr4 [m/s^2]')
    xlabel('t [s]')
    subplot(4,3,10)
    plot(t,ddr8)
    ylabel('ddr8 [m/s^2]')
    xlabel('t [s]')
    subplot(4,3,11)
    plot(t,ddr13)
    ylabel('ddr13 [m/s^2]')
    xlabel('t [s]')

    
% POSITION CONTROL
    F_check2 = zeros(t_size,1);
    f = zeros(t_size,1);
    l = zeros(t_size,1);
    c = zeros(t_size,1);
    for i=1:t_size
        F_check2(i) = (D + r5*exp(j*phi4(i))+r6*exp(j*phi5(i))) - (G + r7*exp(j*(phi6(i)-pi)));
        f(i) = (D + r5*exp(j*phi4(i))+r6*exp(j*phi5(i)));
    end
    I_check2 = zeros(t_size,1);
    for i=1:t_size
        I_check2(i) =(K + r11*exp(j*(phi11(i)-pi)) + r10*exp(j*(phi10(i)+pi))) -(H +r9*exp(j*phi9(i))) ;
        l(i) = (K + r11*exp(j*(phi11(i)-pi)) + r10*exp(j*(phi10(i)+pi)));
    end
    C_check2 = zeros(t_size,1);
    for i=1:t_size
        C_check2(i) =(A + r2*exp(j*phi1(i)) - r3*exp(j*phi2(i))) -(D + r4(i)*exp(j*(phi3(i)+pi)));
        c(i) = (r2 * exp(j*phi1(i)) - r3*exp(j*phi2(i)));
    end
    t1 = t(2:end);      % eerste waarde eruit
    F_check2 = F_check2(2:end,:);
    I_check2 = I_check2(2:end,:);
    C_check2 = C_check2(2:end,:);
    f1 = f(2:end,:);
    l1 = l(2:end,:);
    c1 = c(2:end,:);
    figure('Name','Position control for d, l, q via different paths')
    subplot(321)
    plot(t1,abs(F_check2))
    xlabel('t [s]') 
    ylabel('absolute error d [m]') 
    subplot(323)
    plot(t1,abs(I_check2))
    xlabel('t [s]') 
    ylabel('absolute error l [m]') 
    subplot(325)
    plot(t1,abs(C_check2))
    xlabel('t [s]') 
    ylabel('absolute error c [m]') 
    subplot(322)
    plot(t1,abs(F_check2)./abs(f1))
    xlabel('t [s]') 
    ylabel('relative error d []')
    subplot(324)
    plot(t1,abs(I_check2)./abs(l1))
    xlabel('t [s]') 
    ylabel('relative error l []') 
    subplot(326)
    plot(t1,abs(C_check2)./abs(c1))
    xlabel('t [s]') 
    ylabel('relative error c []') 
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Position control for d, l, q via different paths','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on') 
    
    % VELOCITY CONTROL WITH NUMERICAL DIFFERENTIATION
    dphi5_control = diff(phi5)/Ts;
    dphi11_control = diff(phi11)/Ts;
    dphi3_control = diff(phi3)/Ts;  
    
    
    dphi5_control = [dphi5_control ; dphi5_control(end)];
    dphi11_control = [dphi11_control ; dphi11_control(end)];
    dphi3_control = [dphi3_control ; dphi3_control(end)];  
    
    dphi5_check = dphi5-dphi5_control;
    dphi11_check = dphi11-dphi11_control;
    dphi3_check = dphi3 - dphi3_control;
    
    figure('Name','Velocity control for bar 6, 11, 4 with numerical differentiation')
    subplot(321)
    plot(t,dphi5_check)
    xlabel('t [s]')
    ylabel('absolute error \omega_5 [rad/s]')
    subplot(323)
    plot(t,dphi11_check)
    xlabel('t [s]')
    ylabel('absolute error \omega_11 [rad/s]')
    subplot(325)
    xlabel('t [s]')
    plot(t1,dphi3_check(2:end))
    xlabel('t [s]')
    ylabel('absolute error \omega_3 [rad/s]')
    subplot(322)
    xlabel('t [s]')
    plot(t,dphi11_check./dphi11)
    xlabel('t [s]')
    ylabel('relative error \omega_8 []')   
    subplot(324)
    xlabel('t [s]')
    plot(t,dphi5_check./dphi5)
    xlabel('t [s]')
    ylabel('relative error \omega_1_4 []')   
    subplot(326)
    xlabel('t [s]')
    plot(t,dphi3_check./dphi3)
    xlabel('t [s]')
    ylabel('relative error \omega_5_b []')       
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Velocity control for bar 6, 11, 4 with numerical differentiation','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on')
    
    % ACCELERATION CONTROL WITH NUMERICAL DIFFERENTIATION
    ddphi5_control = diff(dphi5)/Ts;
    ddphi11_control = diff(dphi11)/Ts;
    ddphi3_control = diff(dphi3)/Ts;  
    
    
    ddphi5_control = [ddphi5_control ; ddphi5_control(end)];
    ddphi11_control = [ddphi11_control ; ddphi11_control(end)];
    ddphi3_control = [ddphi3_control ; ddphi3_control(end)];  
    
    ddphi5_check = ddphi5-ddphi5_control;
    ddphi11_check = ddphi11-ddphi11_control;
    ddphi3_check = ddphi3-ddphi3_control;
    
    figure('Name','Acceleration control for bar 6, 11, 4 with numerical differentiation')
    subplot(321)
    plot(t,ddphi5_check)
    xlabel('t [s]')
    ylabel('absolute error \alpha_5 [rad/s^2] ')
    subplot(323)
    plot(t,ddphi11_check)
    xlabel('t [s]')
    ylabel('absolute error \alpha_11 [rad/s^2] ')
    subplot(325)    
    plot(t1,ddphi3_check(2:end))
    xlabel('t [s]')
    ylabel('absolute error \alpha_3[rad/s^2] ')
    subplot(322)
    plot(t,ddphi5_check./ddphi5)
    xlabel('t [s]')
    ylabel('relative error \alpha_5 [] ')
    subplot(324)
    plot(t,ddphi11_check./ddphi11)
    xlabel('t [s]')
    ylabel('relative error \alpha_11 [] ')
    subplot(326)
    plot(t,ddphi3_check./ddphi3)
    xlabel('t [s]')
    ylabel('relative error \alpha_3 [] ')

    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Acceleration control for bar 8, 14, 5 with numerical differentiation','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on')
%     
%   for i=1:t_size  
%   r2_vec=   [r2*cos(phi1)                   r2*sin(phi1)                zeros(size(phi1))];
%   r3_vec=   [r3*cos(phi2)                   r3*sin(phi2)                zeros(size(phi1))];
%   r4_vec=   [r4(i)*cos(phi3)                   r4(i)*sin(phi3)                zeros(size(phi1))];
%   r5_vec=   [r5*cos(phi3 + gamma)           r5*sin(phi3 + gamma)        zeros(size(phi1))];
%   r6_vec=   [r6*cos(2*pi - phi5)                   r6*sin(2*pi -phi5)                zeros(size(phi1))];
%   r7_vec=   [r7*cos(phi6 -pi)                   r7*sin(phi6 -pi)                zeros(size(phi1))];
%   r8_vec=   [r8(i)*cos(phi8)                   r8(i)*sin(phi8)                zeros(size(phi1))];
%   r9_vec=   [r9*cos(phi3 +2*pi - gamma)     r9*sin(phi3 + 2*pi - gamma) zeros(size(phi1))];
%   r10_vec=   [r10*cos(phi10)                 r10*sin(phi10)              zeros(size(phi1))];
%   r11_vec=   [r11*cos(phi11)                 r11*sin(phi11)              zeros(size(phi1))];
%   end
%   %Hoeksnelheden als vectoren
%   omega1 = [zeros(size(phi1))       zeros(size(phi1))   dphi1];
%   omega2 = [zeros(size(phi1))       zeros(size(phi1))   dphi2];
%   omega3 = [zeros(size(phi1))       zeros(size(phi1))   dphi3];
%   omega5 = [zeros(size(phi1))       zeros(size(phi1))   dphi5];
%   omega6 = [zeros(size(phi1))       zeros(size(phi1))   dphi6];
%   omega8 = [zeros(size(phi1))       zeros(size(phi1))   dphi8];
%   omega10 = [zeros(size(phi1))       zeros(size(phi1))   dphi10];
%   omega11= [zeros(size(phi1))       zeros(size(phi1))   dphi11];
%   %Hoekversnellingen als vectoren
%   alpha1 = [zeros(size(phi1))       zeros(size(phi1))   ddphi1];
%   alpha2 = [zeros(size(phi1))       zeros(size(phi1))   ddphi2];
%   alpha3 = [zeros(size(phi1))       zeros(size(phi1))   ddphi3];
%   alpha5 = [zeros(size(phi1))       zeros(size(phi1))   ddphi5];
%   alpha6 = [zeros(size(phi1))       zeros(size(phi1))   ddphi6];
%   alpha8 = [zeros(size(phi1))       zeros(size(phi1))   ddphi8];
%   alpha10 = [zeros(size(phi1))       zeros(size(phi1))   ddphi10];
%   alpha11 = [zeros(size(phi1))       zeros(size(phi1))   ddphi11];
%   
    %Positie Controle
  %De posities van enkele punten worden op twee verschillende manieren
  %uitgerekend en vergeleken
%   
%   %Punt D
%   p_E_1 = r5_vec
%   display r5_vec;
%   p_E_2 = r7_vec-r6_vec
%   display p_E_2;
%   p_E_diff = p_E_1-p_E_2
%   display P_E_diff;
%   p_E_check=sqrt(p_E_diff(:,1).^2+p_E_diff(:,2).^2);
%   
%   %Punt E
%   p_E_1 = r1_vec+r4_vec+r7_vec;
%   p_E_2 = r2_vec+r6_vec+r8_vec;
%   p_E_diff = p_E_1-p_E_2;
%   p_E_check=sqrt(p_E_diff(:,1).^2+p_E_diff(:,2).^2);
%   
%   %Punt F
%   p_F_1 = r1_vec+r3_vec;
%   p_F_2 = r2_vec+r6_vec;
%   p_F_diff = p_F_1-p_F_2;
%   p_F_check=sqrt(p_F_diff(:,1).^2+p_F_diff(:,2).^2);
  
  %Uitzetten van de verschillen tussen de berekende locaties van de punten 
%   figure
%   subplot(1,3,1)
%   plot(t,p_E_check)
%   ylabel('\Delta_p_D [m]')
%   xlabel('t [s]')
%   xlim([0 20*pi])
%   
%   subplot(1,3,2)
%   plot(t,p_E_check)
%   ylabel('\Delta_p_E [m]')
%   xlabel('t [s]')
%   xlim([0 20*pi])
%   
%   subplot(1,3,3)
%   plot(t,p_F_check)
%   ylabel('\Delta_p_F [m]')
%   xlabel('t [s]')
%   xlim([0 20*pi])
%   
%   %SnelheidsControle
%   %Gelijkaardig aan de positiecontrole in concept maar nu worden de
%   %snelheden van punten vergeleken en niet hun posities
%   
%   %Punt D
%   v_D_1=    cross(omega4,r4_vec);
%   v_D_2=    cross(omega2,r2_vec)+cross(omega5,r5_vec);
%   v_D_diff = v_D_1-v_D_2;
%   v_D_check=sqrt(v_D_diff(:,1).^2+v_D_diff(:,2).^2);
%   
%   %Punt E
%   v_E_1=    cross(omega3,r3_vec)+cross(omega8,r8_vec);
%   v_E_2=    cross(omega4,r4_vec)+cross(omega7,r7_vec);
%   v_E_diff = v_E_1-v_E_2;
%   v_E_check=sqrt(v_E_diff(:,1).^2+v_E_diff(:,2).^2);
%   
%   %Punt F
%   v_F_1=    cross(omega3,r3_vec);
%   v_F_2=    cross(omega2,r2_vec)+cross(omega6,r6_vec);
%   v_F_diff = v_F_1-v_F_2;
%   v_F_check=sqrt(v_F_diff(:,1).^2+v_F_diff(:,2).^2);
%   
%   %Uitzetten van de verschillen tussen de berekende snelheden van de punten 
%   figure
%   subplot(1,3,1)
%   plot(t,v_D_check)
%   ylabel('\Delta_v_D [m/s]')
%   xlabel('t [s]')
%   xlim([0 20*pi])
%   
%   subplot(1,3,2)
%   plot(t,v_E_check)
%   ylabel('\Delta_v_E [m/s]')
%   xlabel('t [s]')
%   xlim([0 20*pi])
%   
%   subplot(1,3,3)
%   plot(t,v_F_check)
%   ylabel('\Delta_v_F [m/s]')
%   xlabel('t [s]')
%   xlim([0 20*pi])
%   
%   %Versnellingsanalyse
%   %Gelijkardig als bij de snelheidsanalyse wordt hier de versnelling
%   %gecontroleerd door hem te berekenen op twee verschillende manieren en
%   %het verschil uit te zetten in functie van phi2
%   
%   %Punt D
%   a_D_1=    cross(alpha4,r4_vec)+cross(omega4,cross(omega4,r4_vec));
%   a_D_2=    cross(alpha2,r2_vec)+cross(omega2,cross(omega2,r2_vec))+cross(alpha5,r5_vec)+cross(omega5,cross(omega5,r5_vec));
%   a_D_diff = a_D_1-a_D_2;
%   a_D_check=sqrt(a_D_diff(:,1).^2+a_D_diff(:,2).^2);
%   
%   %Punt E
%   a_E_1=    cross(alpha4,r4_vec)+cross(omega4,cross(omega4,r4_vec))+cross(alpha7,r7_vec)+cross(omega7,cross(omega7,r7_vec));
%   a_E_2=    cross(alpha3,r3_vec)+cross(omega3,cross(omega3,r3_vec))+cross(alpha8,r8_vec)+cross(omega8,cross(omega8,r8_vec));
%   a_E_diff = a_E_1-a_E_2;
%   a_E_check=sqrt(a_E_diff(:,1).^2+a_E_diff(:,2).^2);
%   
%   %Punt F
%   a_F_1=    cross(alpha3,r3_vec)+cross(omega3,cross(omega3,r3_vec));
%   a_F_2=    cross(alpha2,r2_vec)+cross(omega2,cross(omega2,r2_vec))+cross(alpha6,r6_vec)+cross(omega6,cross(omega6,r6_vec));
%   a_F_diff = a_F_1-a_F_2;
%   a_F_check=sqrt(a_F_diff(:,1).^2+a_F_diff(:,2).^2);
%   
%   %Uitzetten van de verschillen tussen de berekende versnellingen van de punten
%   figure
%   subplot(1,3,1)
%   plot(t,a_D_check)
%   ylabel('\Delta_a_D [m/s²]')
%   xlabel('t [s]')
%   xlim([0 20*pi])
%   
%   subplot(1,3,2)
%   plot(t,a_E_check)
%   ylabel('\Delta_a_E [m/s²]')
%   xlabel('t [s]')
%   xlim([0 20*pi])
%   
%   subplot(1,3,3)
%   plot(t,a_F_check)
%   ylabel('\Delta_a_F [m/s²]')
%   xlabel('t [s]')
%   xlim([0 20*pi])
  



 
=======
% 
% 
% % *** create movie ***
% 
% % point P = fixed
% A = 0;
% % point S = fixed
% D = A +0.007746*exp(j*-63.14*pi/180);
% G = A +0.010983*exp(j*-71.4166*pi/180);
% H = A -0.007746*exp(j*63.14*pi/180);
% K = A -0.010983*exp(j*71.4166*pi/180);
% 
% 
% % define which positions we want as frames in our movie
% frames = 200;    % number of frames in movie
% delta = floor(t_size/frames); % time between frames
% index_vec = [1:delta:t_size]';
% 
% % Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% % This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% % axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% % plots.
% x_left = -0.015;
% y_bottom = -0.015;
% x_right = 0.015;
% y_top = 0.005;
% 
% figure(10)
% hold on
% plot([x_left, x_right], [y_bottom, y_top]);
% axis equal;
% movie_axes = axis;   %save current axes into movie_axes
% 
% % draw and save movie frame
% for m=1:length(index_vec)
%     index = index_vec(m);
%     B = A + r2 * exp(j*phi1(index));
%     C = B - r3 * exp(j*phi2(index));
%     
%     E = D + r5*exp(j*phi4(index));
%     F = G + r7*exp(j*(phi6(index)-pi));
% 
%     I = H +r9*exp(j*phi9(index));
%     J = K + r11*exp(j*(phi11(index)-pi));
%     loop = [A B C D E F G F E D C H I J K];
%     figure(10)
%     clf
%     hold on
%     plot(real(loop),imag(loop),'-o')
%     
%     axis(movie_axes);
%     xlabel('x [m]')
%     ylabel('y [m]')
%     set(findobj('type','axes'),'xgrid','on')
%     set(findobj('type','axes'),'ygrid','on')% set axes as in movie_axes
%     Movie(m) = getframe;  % save frame to a variable Film
% end
% 
% % save movie
% save fourbar_movie Movie
% close(10)
% 
% 
% % *** plot figures ***
%     figure('Name','Position: \phi_2 - \r13')
%     subplot(4,3,1)
%     plot(t,phi1)
%     ylabel('\phi_1 [rad]')
%     subplot(4,3,2)
%     plot(t,phi2)
%     ylabel('\phi_2 [rad]')
%     subplot(4,3,3)
%     plot(t,phi3)
%     ylabel('\phi_3 [rad]')
%     subplot(4,3,4)
%     plot(t,phi5)
%     ylabel('\phi_5 [rad]')
%     xlabel('t [s]')
%     subplot(4,3,5)
%     plot(t,phi6)
%     ylabel('\phi_6 [rad]')
%     xlabel('t [s]')
%     subplot(4,3,6)
%     plot(t,phi8)
%     ylabel('\phi_8 [rad]')
%     xlabel('t [s]')
%     subplot(4,3,7)
%     plot(t,phi10)
%     ylabel('\phi_10 [rad]')
%     xlabel('t [s]')
%     subplot(4,3,8)
%     plot(t,phi11)
%     ylabel('\phi_11 [rad]')
%     xlabel('t [s]')
%     subplot(4,3,9)
%     plot(t,r4)
%     ylabel('r4 [m]')
%     xlabel('t [s]')
%     subplot(4,3,10)
%     plot(t,r8)
%     ylabel('r8 [m]')
%     xlabel('t [s]')
%     subplot(4,3,11)
%     plot(t,r13)
%     ylabel('r13 [m]')
%     xlabel('t [s]')
%     
%     figure
%     subplot(4,3,1)
%     plot(t,dphi1)
%     ylabel('d\phi_1 [rad/s]')
%     xlabel('t [s]')
%     subplot(4,3,2)
%     plot(t,dphi2)
%     ylabel('d\phi_2 [rad/s]')
%     xlabel('t [s]')
%     subplot(4,3,3)
%     plot(t,dphi3)
%     ylabel('d\phi_3 [rad/s]')
%     xlabel('t [s]')
%     subplot(4,3,4)
%     plot(t,dphi5)
%     ylabel('d\phi_5 [rad/s]')
%     xlabel('t [s]')
%     subplot(4,3,5)
%     plot(t,dphi6)
%     ylabel('d\phi_6 [rad/s]')
%     xlabel('t [s]')
%     subplot(4,3,6)
%     plot(t,dphi8)
%     ylabel('d\phi_8 [rad/s]')
%     xlabel('t [s]')
%     subplot(4,3,7)
%     plot(t,dphi10)
%     ylabel('d\phi_10 [rad/s]')
%     xlabel('t [s]')
%     subplot(4,3,8)
%     plot(t,dphi11)
%     ylabel('d\phi_11 [rad/s]')
%     xlabel('t [s]')
%     subplot(4,3,9)
%     plot(t,dr4)
%     ylabel('dr4 [m/s]')
%     xlabel('t [s]')
%     subplot(4,3,10)
%     plot(t,dr8)
%     ylabel('dr8 [m/s]')
%     xlabel('t [s]')
%     subplot(4,3,11)
%     plot(t,dr13)
%     ylabel('dr13 [m/s]')
%     xlabel('t [s]')
%     
%     figure
%     subplot(4,3,1)
%     plot(t,ddphi1)
%     ylabel('dd\phi_1 [rad/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,2)
%     plot(t,ddphi2)
%     ylabel('dd\phi_2 [rad/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,3)
%     plot(t,ddphi3)
%     ylabel('dd\phi_3 [rad/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,4)
%     plot(t,ddphi5)
%     ylabel('dd\phi_5 [rad/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,5)
%     plot(t,ddphi6)
%     ylabel('dd\phi_6 [rad/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,6)
%     plot(t,ddphi8)
%     ylabel('dd\phi_8 [rad/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,7)
%     plot(t,ddphi10)
%     ylabel('dd\phi_10 [rad/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,8)
%     plot(t,ddphi11)
%     ylabel('dd\phi_11 [rad/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,9)
%     plot(t,ddr4)
%     ylabel('ddr4 [m/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,10)
%     plot(t,ddr8)
%     ylabel('ddr8 [m/s^2]')
%     xlabel('t [s]')
%     subplot(4,3,11)
%     plot(t,ddr13)
%     ylabel('ddr13 [m/s^2]')
%     xlabel('t [s]')
% 
%     
% % POSITION CONTROL
%     F_check2 = zeros(t_size,1);
%     f = zeros(t_size,1);
%     l = zeros(t_size,1);
%     q = zeros(t_size,1);
%     for i=1:t_size
%         F_check2(i) = (D + r5*exp(j*phi4(i))+r6*exp(j*phi5(i))) - (G + r7*exp(j*(phi6(i)-pi)));
%         f(i) = (D + r5*exp(j*phi4(i))+r6*exp(j*phi5(i)));
%     end
%     L_check2 = zeros(t_size,1);
%     for i=1:t_size
%         L_check2(i) =(K + r11*exp(j*(phi11(i)-pi)) + r10*exp(j*(phi10(i)+pi))) -(H +r9*exp(j*phi9(i))) ;
%         l(i) = (K + r11*exp(j*(phi11(i)-pi)) + r10*exp(j*(phi10(i)+pi)));
%     end
%     Q_check2 = zeros(t_size,1);
%     for i=1:t_size
%         Q_check2(i) =(A + r2*exp(j*phi1(i)) - r3*exp(j*phi2(i))) -(D + r4(i)*exp(j*(phi3(i)+pi)));
%         q(i) = (r2 * exp(j*phi1(i)) - r3*exp(j*phi2(i)));
%     end
%     t1 = t(2:end);      % eerste waarde eruit
%     F_check2 = F_check2(2:end,:);
%     L_check2 = L_check2(2:end,:);
%     Q_check2 = Q_check2(2:end,:);
%     f1 = f(2:end,:);
%     l1 = l(2:end,:);
%     q1 = q(2:end,:);
%     figure('Name','Position control for d, l, q via different paths')
%     subplot(321)
%     %plot(t,abs(D_check))
%     plot(t1,abs(F_check2))
%     xlabel('t [s]') 
%     ylabel('absolute error d [m]') 
%     subplot(323)
%     %plot(t,abs(L_check))
%     plot(t1,abs(L_check2))
%     xlabel('t [s]') 
%     ylabel('absolute error l [m]') 
%     subplot(325)
%     plot(t1,abs(Q_check2))
%     xlabel('t [s]') 
%     ylabel('absolute error c [m]') 
%     subplot(322)
%     plot(t1,abs(F_check2)./abs(f1))
%     xlabel('t [s]') 
%     ylabel('relative error d []')
%     subplot(324)
%     %plot(t,abs(L_check))
%     plot(t1,abs(L_check2)./abs(l1))
%     xlabel('t [s]') 
%     ylabel('relative error l []') 
%     subplot(326)
%     %plot(t,abs(L_check))
%     plot(t1,abs(Q_check2)./abs(q1))
%     xlabel('t [s]') 
%     ylabel('relative error c []') 
%     
%     ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
%     1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%     text(0.5, 1,'\bf Position control for d, l, q via different paths','HorizontalAlignment','center','VerticalAlignment', 'top')
%     set(findobj('type','axes'),'xgrid','on')
%     set(findobj('type','axes'),'ygrid','on') 
%     
%     % VELOCITY CONTROL WITH NUMERICAL DIFFERENTIATION
%     dphi5_control = diff(phi5)/Ts;
%     dphi11_control = diff(phi11)/Ts;
%     dphi3_control = diff(phi3)/Ts;  
%     
%     
%     dphi5_control = [dphi5_control ; dphi5_control(end)];
%     dphi11_control = [dphi11_control ; dphi11_control(end)];
%     dphi3_control = [dphi3_control ; dphi3_control(end)];  
%     
%     dphi5_check = dphi5-dphi5_control;
%     dphi11_check = dphi11-dphi11_control;
%     dphi3_check = dphi3 - dphi3_control;
%     
%     figure('Name','Velocity control for bar 8, 14, 5 with numerical differentiation')
%     subplot(321)
%     plot(t,dphi5_check)
%     xlabel('t [s]')
%     ylabel('absolute error \omega_8 [rad/s]')
%     subplot(323)
%     plot(t,dphi11_check)
%     xlabel('t [s]')
%     ylabel('absolute error \omega_1_4 [rad/s]')
%     subplot(325)
%     xlabel('t [s]')
%     plot(t1,dphi3_check(2:end))
%     xlabel('t [s]')
%     ylabel('absolute error \omega_5_b [rad/s]')
%     subplot(322)
%     xlabel('t [s]')
%     plot(t,dphi11_check./dphi11)
%     xlabel('t [s]')
%     ylabel('relative error \omega_8 []')   
%     subplot(324)
%     xlabel('t [s]')
%     plot(t,dphi5_check./dphi5)
%     xlabel('t [s]')
%     ylabel('relative error \omega_1_4 []')   
%     subplot(326)
%     xlabel('t [s]')
%     plot(t,dphi3_check./dphi3)
%     xlabel('t [s]')
%     ylabel('relative error \omega_5_b []')       
%     ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
%     1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%     text(0.5, 1,'\bf Velocity control for bar 8, 14, 5 with numerical differentiation','HorizontalAlignment','center','VerticalAlignment', 'top')
%     set(findobj('type','axes'),'xgrid','on')
%     set(findobj('type','axes'),'ygrid','on')
%     
%     % ACCELERATION CONTROL WITH NUMERICAL DIFFERENTIATION
%     ddphi5_control = diff(dphi5)/Ts;
%     ddphi11_control = diff(dphi11)/Ts;
%     ddphi3_control = diff(dphi3)/Ts;  
%     
%     
%     ddphi5_control = [ddphi5_control ; ddphi5_control(end)];
%     ddphi11_control = [ddphi11_control ; ddphi11_control(end)];
%     ddphi3_control = [ddphi3_control ; ddphi3_control(end)];  
%     
%     ddphi5_check = ddphi5-ddphi5_control;
%     ddphi11_check = ddphi11-ddphi11_control;
%     ddphi3_check = ddphi3 - ddphi3_control;
%     
%     figure('Name','Acceleration control for bar 8, 14, 5 with numerical differentiation')
%     subplot(321)
%     plot(t,ddphi5_check)
%     xlabel('t [s]')
%     ylabel('absolute error \alpha_8 [rad/s^2] ')
%     subplot(323)
%     plot(t,ddphi11_check)
%     xlabel('t [s]')
%     ylabel('absolute error \alpha_1_4 [rad/s^2] ')
%     subplot(325)    
%     plot(t1,ddphi3_check(2:end))
%     xlabel('t [s]')
%     ylabel('absolute error \alpha_5_b [rad/s^2] ')
%     subplot(322)
%     plot(t,ddphi5_check./ddphi5)
%     xlabel('t [s]')
%     ylabel('relative error \alpha_8 [] ')
%     subplot(324)
%     plot(t,ddphi11_check./ddphi11)
%     xlabel('t [s]')
%     ylabel('relative error \alpha_1_4 [] ')
%     subplot(326)
%     plot(t,ddphi3_check./ddphi3)
%     xlabel('t [s]')
%     ylabel('relative error \alpha_5_b [] ')
% 
%     ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
%     1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%     text(0.5, 1,'\bf Acceleration control for bar 8, 14, 5 with numerical differentiation','HorizontalAlignment','center','VerticalAlignment', 'top')
%     set(findobj('type','axes'),'xgrid','on')
%     set(findobj('type','axes'),'ygrid','on')
% 
% 
% 
%  
>>>>>>> 5a367f48cfdb2ccf3003d7edfa2eafda37131708
