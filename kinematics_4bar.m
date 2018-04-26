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

function [phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,dphi3,dphi4,dphi6,dphi7,dphi8,dphi10,dphi11,ddphi3,ddphi4,ddphi6,ddphi7,ddphi8,ddphi10,ddphi11,r8,dr8,ddr8,r13,dr13,ddr13,r4,dr4,ddr4] = kinematics_4bar(r2,r3,r5,r6,r7,r9,r10,r11,r14x,r14y,r47y,r18x,r18y,r811y,phi2,dphi2,ddphi2,phi3_init,phi4_init,phi6_init,phi7_init,phi8_init,phi10_init,phi11_init,r13_init,r4_init,r8_init,t,fig_kin_4bar,fig_kin_check);

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
phi3 = zeros(size(t));
phi4 = zeros(size(t));
phi6 = zeros(size(t));
phi7 = zeros(size(t));
phi8 = zeros(size(t));
phi10 = zeros(size(t));
phi11 = zeros(size(t));
dphi3 = zeros(size(t));
dphi4 = zeros(size(t));
dphi6 = zeros(size(t));
dphi7 = zeros(size(t));
dphi8 = zeros(size(t));
dphi10 = zeros(size(t));
dphi11 = zeros(size(t));
ddphi3 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi7 = zeros(size(t));
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

phi5 = zeros(size(t));
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
    [x, fval, exitflag]=fsolve('loop_closure_eqs',[phi3_init phi4_init phi6_init phi7_init phi8_init phi10_init phi11_init r13_init r4_init r8_init]',optim_options,phi2(k),r2,r3,r5,r6,r7,r9,r10,r11,gamma,r14x,r14y,r47y,r18x,r18y,r811y);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % save results of fsolve
    phi3(k)=x(1);
    phi4(k)=x(2);
    phi6(k)=x(3);
    phi7(k)=x(4);
    phi8(k)=x(5);
    phi10(k)=x(6);
    phi11(k)=x(7);
    r13(k)=x(8);
    r4(k)=x(9);
    r8(k)=x(10);
    
    phi5(k) = phi4(k) + gamma;
    phi9(k) = phi8(k) + 2*pi - gamma;
    
    % *** velocity analysis ***
    
    A = [r3*sin(phi3(k)),0,0,0,0,0,0,0,0,0;
        -r3*cos(phi3(k)),0,0,0,0,0,0,1,0,0;
        0,r4(k)*sin(phi4(k)),0,0,0,0,0,0,-cos(phi4(k)),0;
        0,-r4(k)*cos(phi4(k)),0,0,0,0,0,1,-sin(phi4(k)),0;
        0,0,0,0,r8(k)*sin(phi8(k)),0,0,0,0,-cos(phi8(k));
        0,0,0,0,-r8(k)*cos(phi8(k)),0,0,1,0,-sin(phi8(k));
        0,-r5*sin(phi5(k)),-r6*sin(phi6(k)),-r7*sin(phi7(k)),0,0,0,0,0,0;
        0,r5*cos(phi5(k)),r6*cos(phi6(k)),r7*cos(phi7(k)),0,0,0,0,0,0;
        0,0,0,0,-r9*sin(phi9(k)),-r10*sin(phi10(k)),-r11*sin(phi11(k)),0,0,0;
        0,0,0,0,r9*cos(phi9(k)),r10*cos(phi10(k)),r11*cos(phi11(k)),0,0,0];
        
    B = [+r2*sin(phi2(k))*dphi2(k);
        -r2*cos(phi2(k))*dphi2(k);
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
    dphi3(k)=x(1);
    dphi4(k)=x(2);
    dphi6(k)=x(3);
    dphi7(k)=x(4);
    dphi8(k)=x(5);
    dphi10(k)=x(6);
    dphi11(k)=x(7);
    dr13(k)=x(8);
    dr4(k)=x(9);
    dr8(k)=x(10);
    
    
    % *** acceleration analysis ***
    
        
      A = [r3*sin(phi3(k)),0,0,0,0,0,0,0,0,0;
          -r3*cos(phi3(k)),0,0,0,0,0,0,1,0,0;
          0,r4(k)*sin(phi4(k)),0,0,0,0,0,0,-cos(phi4(k)),0;
          0,-r4(k)*cos(phi4(k)),0,0,0,0,0,1,-sin(phi4(k)),0;
          0,0,0,0,r8(k)*sin(phi8(k)),0,0,0,0,-cos(phi8(k));
          0,0,0,0,-r8(k)*cos(phi8(k)),0,0,1,0,-sin(phi8(k));
          0,-r5*sin(phi5(k)),-r6*sin(phi6(k)),-r7*sin(phi7(k)),0,0,0,0,0,0;
          0,r5*cos(phi5(k)),r6*cos(phi6(k)),r7*cos(phi7(k)),0,0,0,0,0,0;
          0,0,0,0,-r9*sin(phi9(k)),-r10*sin(phi10(k)),-r11*sin(phi11(k)),0,0,0;
          0,0,0,0,r9*cos(phi9(k)),r10*cos(phi10(k)),r11*cos(phi11(k)),0,0,0];
    
      B = [r2*cos(phi2(k))*dphi2(k)^2+r2*sin(phi2(k))*ddphi2(k)-r3*cos(phi3(k))*dphi3(k)^2;
          -r2*cos(phi2(k))*ddphi2(k)+r2*sin(phi2(k))*dphi2(k)^2-r3*sin(phi3(k))*dphi3(k)^2;
          -dphi4(k)^2*r4(k)*cos(phi4(k))-dphi4(k)*sin(phi4(k))*dr4(k)-dr4(k)*dphi4(k)*sin(phi4(k));
          -dphi4(k)^2*r4(k)*sin(phi4(k))+dphi4(k)*cos(phi4(k))*dr4(k)+dr4(k)*dphi4(k)*cos(phi4(k));
          -dphi8(k)^2*r8(k)*cos(phi8(k))-dphi8(k)*sin(phi8(k))*dr8(k)-dr8(k)*dphi8(k)*sin(phi8(k));
          -dphi8(k)^2*r8(k)*sin(phi8(k))+dphi8(k)*cos(phi8(k))*dr8(k)+dr8(k)*dphi8(k)*cos(phi8(k));
          r5*cos(phi5(k))*dphi4(k)^2+r6*cos(phi6(k))*dphi6(k)^2+r7*cos(phi7(k))*dphi7(k)^2;
          r5*sin(phi5(k))*dphi4(k)^2+r6*sin(phi6(k))*dphi6(k)^2+r7*sin(phi7(k))*dphi7(k)^2;
          r9*cos(phi9(k))*dphi8(k)^2+r10*cos(phi10(k))*dphi10(k)^2+r11*cos(phi11(k))*dphi11(k)^2;
          r9*sin(phi9(k))*dphi8(k)^2+r10*sin(phi10(k))*dphi10(k)^2+r11*sin(phi11(k))*dphi11(k)^2
          ];
    
    x = A\B;
    % save results
    ddphi3(k)=x(1);
    ddphi4(k)=x(2);
    ddphi6(k)=x(3);
    ddphi7(k)=x(4);
    ddphi8(k)=x(5);
    ddphi10(k)=x(6);
    ddphi11(k)=x(7);
    ddr13(k)=x(8);
    ddr4(k)=x(9);
    ddr8(k)=x(10);
    
    
    % *** calculate initial values for next iteration step ***
    phi3_init = phi3(k)+Ts*dphi3(k);
    phi4_init = phi4(k)+Ts*dphi4(k);
    phi6_init = phi6(k)+Ts*dphi6(k);
    phi7_init = phi7(k)+Ts*dphi7(k);
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
if fig_kin_4bar
% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    B = A + r2 * exp(j*phi2(index));
    C = B - r3 * exp(j*phi3(index));
    
    E = D + r5*exp(j*phi5(index));
    F = G + r7*exp(j*(phi7(index)-pi));

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
save fourbar_movie Movie;
close(10);


    
    %plot assembly at a certain timestep 
    B = A + r2 * exp(j*phi2(1));
    C = B - r3 * exp(j*phi3(1));
    
    E = D + r5*exp(j*phi5(1));
    F = G + r7*exp(j*(phi7(1)-pi));

    I = H +r9*exp(j*phi9(1));
    J = K + r11*exp(j*(phi11(1)-pi));
    
    figure
    assembly=[A B C D E F G F E D C H I J K];
    plot(real(assembly),imag(assembly),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal


dphi5 = dphi4;
ddphi5 = ddphi4;
dphi9 = dphi8;
ddphi9 = ddphi8;
% *** plot figures ***
    figure('Name','Position: \phi_2 - \phi_10')
    subplot(3,3,1)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(3,3,2)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(3,3,3)
    plot(t,phi4)
    ylabel('\phi_4 [rad]')
    subplot(3,3,4)
    plot(t,phi5)
    ylabel('\phi_5 [rad]')
    xlabel('t [s]')
    subplot(3,3,5)
    plot(t,phi6)
    ylabel('\phi_6 [rad]')
    xlabel('t [s]')
    subplot(3,3,6)
    plot(t,phi7)
    ylabel('\phi_7 [rad]')
    xlabel('t [s]')
    subplot(3,3,7)
    plot(t,phi8)
    ylabel('\phi_8 [rad]')
    xlabel('t [s]')
    subplot(3,3,8)
    plot(t,phi9)
    ylabel('\phi_9 [rad]')
    xlabel('t [s]')
    subplot(3,3,9)
    plot(t,phi10)
    ylabel('\phi_10 [rad]')
    xlabel('t [s]')
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'Position: \phi_2 - \phi_10','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on') 
    
    figure('Name','Position: \phi_11 - r13')
    subplot(4,1,1)
    plot(t,phi11)
    ylabel('\phi_11 [rad]')
    xlabel('t [s]')
    subplot(4,1,2)
    plot(t,r4)
    ylabel('r4 [m]')
    xlabel('t [s]')
    subplot(4,1,3)
    plot(t,r8)
    ylabel('r8 [m]')
    xlabel('t [s]')
    subplot(4,1,4)
    plot(t,r13)
    ylabel('r13 [m]')
    xlabel('t [s]')
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'Position: \phi_11 - r13','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on') 
    
    figure('Name','Velocity: d\phi_2 - d\phi_10')
    subplot(3,3,1)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    xlabel('t [s]')
    subplot(3,3,2)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    xlabel('t [s]')
    subplot(3,3,3)
    plot(t,dphi4)
    ylabel('d\phi_4 [rad/s]')
    xlabel('t [s]')
    subplot(3,3,4)
    plot(t,dphi5)
    ylabel('d\phi_5 [rad/s]')
    xlabel('t [s]')
    subplot(3,3,5)
    plot(t,dphi6)
    ylabel('d\phi_6 [rad/s]')
    xlabel('t [s]')
    subplot(3,3,6)
    plot(t,dphi7)
    ylabel('d\phi_7 [rad/s]')
    xlabel('t [s]')
    subplot(3,3,7)
    plot(t,dphi8)
    ylabel('d\phi_8 [rad/s]')
    xlabel('t [s]')
    subplot(3,3,8)
    plot(t,dphi9)
    ylabel('d\phi_9 [rad/s]')
    xlabel('t [s]')
    subplot(3,3,9)
    plot(t,dphi10)
    ylabel('d\phi_10 [rad/s]')
    xlabel('t [s]')
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'Velocity: d\phi_2 - d\phi_10','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on') 
    
    figure('Name','Velocity: d\phi_11 - dr13')
    subplot(4,1,1)
    plot(t,dphi11)
    ylabel('d\phi_11 [rad/s]')
    xlabel('t [s]')
    subplot(4,1,2)
    plot(t,dr4)
    ylabel('dr4 [m/s]')
    xlabel('t [s]')
    subplot(4,1,3)
    plot(t,dr8)
    ylabel('dr8 [m/s]')
    xlabel('t [s]')
    subplot(4,1,4)
    plot(t,dr13)
    ylabel('dr13 [m/s]')
    xlabel('t [s]')
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'Velocity: d\phi_11 - dr13','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on') 
    
    figure('Name','Acceleration: dd\phi_2 - dd\phi_10')
    subplot(3,3,1)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    xlabel('t [s]')
    subplot(3,3,2)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    xlabel('t [s]')
    subplot(3,3,3)
    plot(t,ddphi4)
    ylabel('dd\phi_4 [rad/s^2]')
    xlabel('t [s]')
    subplot(3,3,4)
    plot(t,ddphi5)
    ylabel('dd\phi_5 [rad/s^2]')
    xlabel('t [s]')
    subplot(3,3,5)
    plot(t,ddphi6)
    ylabel('dd\phi_6 [rad/s^2]')
    xlabel('t [s]')
    subplot(3,3,6)
    plot(t,ddphi7)
    ylabel('dd\phi_7 [rad/s^2]')
    xlabel('t [s]')
    subplot(3,3,7)
    plot(t,ddphi8)
    ylabel('dd\phi_8 [rad/s^2]')
    xlabel('t [s]')
    subplot(3,3,8)
    plot(t,ddphi9)
    ylabel('dd\phi_9 [rad/s^2]')
    xlabel('t [s]')
    subplot(3,3,9)
    plot(t,ddphi10)
    ylabel('dd\phi_10 [rad/s^2]')
    xlabel('t [s]')
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'Acceleration: dd\phi_2 - dd\phi_10','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on') 
    
    figure('Name','Acceleration: dd\phi_11 - ddr13')
    subplot(4,1,1)
    plot(t,ddphi11)
    ylabel('dd\phi_11 [rad/s^2]')
    xlabel('t [s]')
    subplot(4,1,2)
    plot(t,ddr4)
    ylabel('ddr4 [m/s^2]')
    xlabel('t [s]')
    subplot(4,1,3)
    plot(t,ddr8)
    ylabel('ddr8 [m/s^2]')
    xlabel('t [s]')
    subplot(4,1,4)
    plot(t,ddr13)
    ylabel('ddr13 [m/s^2]')
    xlabel('t [s]')
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'Acceleration: dd\phi_11 - ddr13','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on') 
end

if fig_kin_check
    
% POSITION CONTROL STARTING FROM DIFFERENT POINTS
    F_check = zeros(t_size,1);
    f = zeros(t_size,1);
    l = zeros(t_size,1);
    c = zeros(t_size,1);
    for i=1:t_size
        F_check(i) = (D + r5*exp(j*phi5(i))+r6*exp(j*phi6(i))) - (G + r7*exp(j*(phi7(i)-pi)));
        f(i) = (D + r5*exp(j*phi5(i))+r6*exp(j*phi6(i)));
    end
    I_check = zeros(t_size,1);
    for i=1:t_size
        I_check(i) =(K + r11*exp(j*(phi11(i)-pi)) + r10*exp(j*(phi10(i)+pi))) -(H +r9*exp(j*phi9(i))) ;
        l(i) = (K + r11*exp(j*(phi11(i)-pi)) + r10*exp(j*(phi10(i)+pi)));
    end
    C_check = zeros(t_size,1);
    for i=1:t_size
        C_check(i) =(A + r2*exp(j*phi2(i)) - r3*exp(j*phi3(i))) -(D + r4(i)*exp(j*(phi4(i)+pi)));
        c(i) = (r2 * exp(j*phi2(i)) - r3*exp(j*phi3(i)));
    end
    t1 = t(2:end);      % eerste waarde eruit
    F_check = F_check(2:end,:);
    I_check = I_check(2:end,:);
    C_check = C_check(2:end,:);
    f1 = f(2:end,:);
    l1 = l(2:end,:);
    c1 = c(2:end,:);
    figure('Name','Position control for F, I, C starting from two different points')
    subplot(321)
    plot(t1,abs(F_check))
    xlabel('t [s]') 
    ylabel('absolute error f [m]') 
    subplot(323)
    plot(t1,abs(I_check))
    xlabel('t [s]') 
    ylabel('absolute error i [m]') 
    subplot(325)
    plot(t1,abs(C_check))
    xlabel('t [s]') 
    ylabel('absolute error c [m]') 
    subplot(322)
    plot(t1,abs(F_check)./abs(f1))
    xlabel('t [s]') 
    ylabel('relative error f []')
    subplot(324)
    plot(t1,abs(I_check)./abs(l1))
    xlabel('t [s]') 
    ylabel('relative error i []') 
    subplot(326)
    plot(t1,abs(C_check)./abs(c1))
    xlabel('t [s]') 
    ylabel('relative error c []') 
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Position control for F, I, C starting from two different points','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on') 
    
    % VELOCITY CONTROL WITH NUMERICAL DIFFERENTIATION
    dphi6_control = diff(phi6)/Ts;
    dphi10_control = diff(phi10)/Ts;
    dphi3_control = diff(phi4)/Ts;  
    
    
    dphi6_control = [dphi6_control ; dphi6_control(end)];
    dphi10_control = [dphi10_control ; dphi10_control(end)];
    dphi3_control = [dphi3_control ; dphi3_control(end)];  
    
    dphi6_check = dphi6-dphi6_control;
    dphi10_check = dphi10-dphi10_control;
    dphi3_check = dphi3 - dphi3_control;
    
    figure('Name','Velocity control for phi 6, 10 and 3 with numerical differentiation')
    subplot(321)
    plot(t,dphi6_check)
    xlabel('t [s]')
    ylabel('absolute error \omega_6 [rad/s]')
    subplot(323)
    plot(t,dphi10_check)
    xlabel('t [s]')
    ylabel('absolute error \omega_10 [rad/s]')
    subplot(325)
    xlabel('t [s]')
    plot(t,dphi3_check)
    xlabel('t [s]')
    ylabel('absolute error \omega_3 [rad/s]')
    subplot(322)
    xlabel('t [s]')
    plot(t,dphi10_check./dphi10)
    xlabel('t [s]')
    ylabel('relative error \omega_10 []')   
    subplot(324)
    xlabel('t [s]')
    plot(t,dphi6_check./dphi6)
    xlabel('t [s]')
    ylabel('relative error \omega_6 []')   
    subplot(326)
    xlabel('t [s]')
    plot(t,dphi3_check./dphi3)
    xlabel('t [s]')
    ylabel('relative error \omega_3 []')       
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Velocity control for phi 6, 10 and 3 with numerical differentiation','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on')
    
    % ACCELERATION CONTROL WITH NUMERICAL DIFFERENTIATION
    ddphi6_control = diff(dphi6)/Ts;
    ddphi10_control = diff(dphi10)/Ts;
    ddphi3_control = diff(dphi3)/Ts;  
    
    
    ddphi6_control = [ddphi6_control ; ddphi6_control(end)];
    ddphi10_control = [ddphi10_control ; ddphi10_control(end)];
    ddphi3_control = [ddphi3_control ; ddphi3_control(end)];  
    
    ddphi6_check = ddphi6-ddphi6_control;
    ddphi10_check = ddphi10-ddphi10_control;
    ddphi3_check = ddphi3-ddphi3_control;
    
    figure('Name','Acceleration control for phi 6, 10 and 3 with numerical differentiation')
    subplot(321)
    plot(t,ddphi6_check)
    xlabel('t [s]')
    ylabel('absolute error \alpha_6 [rad/s^2] ')
    subplot(323)
    plot(t,ddphi10_check)
    xlabel('t [s]')
    ylabel('absolute error \alpha_10 [rad/s^2] ')
    subplot(325)    
    plot(t,ddphi3_check)
    xlabel('t [s]')
    ylabel('absolute error \alpha_3[rad/s^2] ')
    subplot(322)
    plot(t,ddphi6_check./ddphi6)
    xlabel('t [s]')
    ylabel('relative error \alpha_6 [] ')
    subplot(324)
    plot(t,ddphi10_check./ddphi10)
    xlabel('t [s]')
    ylabel('relative error \alpha_10 [] ')
    subplot(326)
    plot(t,ddphi3_check./ddphi3)
    xlabel('t [s]')
    ylabel('relative error \alpha_3 [] ')

    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Acceleration control for phi 6, 10, 3 with numerical differentiation','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on')
    
    % VELOCITY CONTROL USING MULTIPLE POINTS
    % F
    dphi5=dphi4;
    ddphi5=ddphi4;
    dphi9=dphi8;
    ddphi9=ddphi8;
    D_E_vect = [r5*cos(phi5) r5*sin(phi5) zeros(size(phi2))];
    E_F_vect = [r6*cos(phi6) r6*sin(phi6) zeros(size(phi2))];
    omega5 = [zeros(size(phi2)) zeros(size(phi2)) dphi5];
    omega6 = [zeros(size(phi2)) zeros(size(phi2)) dphi6];
    omega7 = [zeros(size(phi2)) zeros(size(phi2)) dphi7];
    G_F_vect = [r7*cos(phi7-pi) r7*sin(phi7-pi) zeros(size(phi2))];
    vF = cross(omega5,D_E_vect) + cross(omega6,E_F_vect);
    vF_check = cross(omega7,G_F_vect);
    vF_diff = zeros(length(vF),1);
    % J    
    H_I_vect = [r9*cos(phi9) r9*sin(phi9) zeros(size(phi2))];
    I_J_vect = [r10*cos(phi10) r10*sin(phi10) zeros(size(phi2))];
    omega9 = [zeros(size(phi2)) zeros(size(phi2)) dphi9];
    omega10 = [zeros(size(phi2)) zeros(size(phi2)) dphi10];
    omega11 = [zeros(size(phi2)) zeros(size(phi2)) dphi11];
    K_J_vect = [r11*cos(phi11-pi) r11*sin(phi11-pi) zeros(size(phi2))];
    vJ = cross(omega9,H_I_vect) + cross(omega10,I_J_vect);
    vJ_check = cross(omega11,K_J_vect);
    vJ_diff = zeros(length(vJ),1);

    for i=1:length(vF)
        vF_diff(i) = norm(vF(i,:)) - norm(vF_check(i,:));
    end
    for i=1:length(vJ)
        vJ_diff(i) = norm(vJ(i,:)) - norm(vJ_check(i,:));
    end
    figure('Name','Velocity control for joints F and J via different paths')
    subplot(221)
    plot(t,vF_diff);
    xlabel('t [s]')
    ylabel('absolute error v_f [m/s] ')
    subplot(222)
    plot(t,vF_diff./norm(vF));
    xlabel('t [s]')
    ylabel('relative error v_f [] ')
    subplot(223)
    plot(t,vJ_diff);
    xlabel('t [s]')
    ylabel('absolute error v_J [m/s] ')
    subplot(224)
    plot(t,vJ_diff./norm(vJ));
    xlabel('t [s]')
    ylabel('relative error v_J [] ')

    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'Velocity control for joints F and J via different paths','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on')
    % ACCELERATION CONTROL USING MULTIPLE POINTS
    % F    
    alpha6 = [zeros(size(phi2)) zeros(size(phi2)) ddphi6];
    alpha7 = [zeros(size(phi2)) zeros(size(phi2)) ddphi7];
    alpha5 = [zeros(size(phi2)) zeros(size(phi2)) ddphi5];
    aF = cross(omega5,cross(omega5,D_E_vect)) + cross(omega6,cross(omega6,E_F_vect)) + cross(alpha5,D_E_vect) + cross(alpha6,E_F_vect);
    aF_check = cross(omega7,cross(omega7,G_F_vect)) + cross(alpha7,G_F_vect);
    aF_diff = zeros(length(aF),1);
    % J
    alpha9 = [zeros(size(phi2)) zeros(size(phi2)) ddphi9];
    alpha10 = [zeros(size(phi2)) zeros(size(phi2)) ddphi10];
    alpha11 = [zeros(size(phi2)) zeros(size(phi2)) ddphi11];
    aJ = cross(omega9,cross(omega9,H_I_vect)) + cross(omega10,cross(omega10,I_J_vect)) + cross(alpha9,H_I_vect) + cross(alpha10,I_J_vect);
    aJ_check = cross(omega11,cross(omega11,K_J_vect)) + cross(alpha11,K_J_vect);
    aJ_diff = zeros(length(aJ),1);  
    for i=1:length(aF)
        aF_diff(i) = norm(aF(i,:)) - norm(aF_check(i,:));
    end
    for i=1:length(aJ)
        aJ_diff(i) = norm(aJ(i,:)) - norm(aJ_check(i,:));
    end
    figure('Name','Acceleration control for joints F and J via different paths')
    subplot(221)
    plot(t,aF_diff);
    xlabel('t [s]')
    ylabel('absolute error a_f [m/s] ')
    subplot(222)
    plot(t,aF_diff./norm(aF));
    xlabel('t [s]')
    ylabel('relative error a_f [] ')
    subplot(223)
    plot(t,aJ_diff);
    xlabel('t [s]')
    ylabel('absolute error a_J [m/s] ')
    subplot(224)
    plot(t,aJ_diff./norm(aJ));
    xlabel('t [s]')
    ylabel('relative error a_f [] ')
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'Acceleration control for joints F and J via different paths','HorizontalAlignment','center','VerticalAlignment', 'top')
    set(findobj('type','axes'),'xgrid','on')
    set(findobj('type','axes'),'ygrid','on')
end
    
