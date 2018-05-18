function eps = Dynamics(zeta,ks,theta1,Omega_rad,Mass,S,normalForce,Pressure_Angle)

disp("start Single rise analysis")
%Rise from 200 till 280 degrees, dwell from 280 till 60 degrees
rise = 200;
dwell = 280;
d_end = 60;
t1 = ((dwell - rise)*pi)/(180*Omega_rad);
tau_end = (d_end-rise+360)/(dwell-rise);
lambda = 0.75/zeta;

kf = Mass*((2*pi*lambda)/(t1))^2 - ks*10^3;

tau = 0:(tau_end/22000):tau_end;

Rise = zeros(size(tau));
Rise(1:16001) = S(20000:36000);
Rise(16002:22001) = S(1:6000);
theta = (Rise)./30;

% The system doesn't start from lift and speed equal to zero
numerator = (2*pi*lambda)^2;
denominator = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
theta0 = 1;
theta_dot0 = 0;
[A,B,C,D] = tf2ss(numerator,denominator);
X0 = [1/C(2)*theta_dot0; 1/C(2)*theta0];

gamma_num = lsim(A,B,C,D, theta, tau, X0);

figure
subplot(2,2,1)
plot(tau,theta)
xlabel('tau [-]')
ylabel('theta [-]')
subplot(2,2,2)
plot(tau,gamma_num)
xlabel('tau [-]')
ylabel('gamma_num [-]')
subplot(2,2,3)
plot(tau, theta-gamma_num.')
xlabel('tau [-]')
ylabel('theta - gamma [-]')

%Calculating the amplitudes
x_0 = gamma_num(8001);
lambda_d = lambda * sqrt(1 - zeta^2);
dx_0 = ((gamma_num(8002) - gamma_num(8000))/(2*0.000125));
A_num = sqrt(((x_0*2*pi*lambda_d)^2 + (dx_0 + zeta*2*pi*lambda*x_0)^2)/((2*pi*lambda_d)^2));
A_approx = (((2*pi)^2)/((2*pi*lambda)^3))*sqrt(1/(1-zeta^2));

eps = (A_num-A_approx)/A_num;

disp("start Multi rise analysis")
%Multi
wn = sqrt((ks*10^3 + kf)/Mass);
tn = (2*pi)/wn;
lambda_tilde = 2/tn;
numerator2 = (2*pi*lambda_tilde)^2;
denominator2 = [1, 2*zeta*(2*pi*lambda_tilde), (2*pi*lambda_tilde)^2];
sys2 = tf(numerator2, denominator2);

% construct theta
tau_multi = theta1/(2*pi);
theta_multi1 = 0.015*((tau_multi - 60/360)/(60/360) - sin(pi*(tau_multi - 60/360)/(60/360))/(pi));
theta_multi2 = 0.015 + 0.015*((tau_multi - 120/360)/(60/360) + sin(pi*(tau_multi - 120/360)/(60/360))/(pi));
theta_multi3 = 0.03*(1 - (tau_multi - 200/360)/(80/360) + sin(2*pi*(tau_multi - 200/360)/(80/360))/(2*pi));

theta_multi = zeros(size(tau_multi));
theta_multi(6001:12001) = theta_multi1(6001:12001);
theta_multi(12001:18001) = theta_multi2(12001:18001);
theta_multi(18001:20001) = 0.03;
theta_multi(20001:28001) = theta_multi3(20001:28001);

theta_multi = theta_multi/0.03;
gamma_num_multi = lsim(sys2, theta_multi, tau_multi); %system starts from rest

figure
plot(tau_multi,theta_multi)
xlabel('tau_multi [-]')
ylabel('theta_multi [-]')

% Fourier series
N = 100;
[a, b] = Fseries(tau_multi, theta_multi, N);

% compute analytical result with fourier coefficients
c = zeros(N, 1);
d = zeros(N, 1);

for k = 1:N
    c(k) = (-2*zeta*k*b(k)/lambda_tilde + a(k+1)*(1 - ((k/lambda_tilde)^2)))/(((2*zeta*k/lambda_tilde)^2) + ((1 - ((k/lambda_tilde)^2))^2));
    d(k) = (2*zeta*k*a(k+1)/lambda_tilde + b(k)*(1 - ((k/lambda_tilde)^2)))/(((2*zeta*k/lambda_tilde)^2) + ((1 - ((k/lambda_tilde)^2))^2));
end

gamma_anal_multi = 0.5*a(1)*ones(size(tau_multi));

for k = 1:N
    gamma_anal_multi = gamma_anal_multi + c(k)*cos(2*pi*k*(tau_multi+0.5)) + d(k)*sin(2*pi*k*(tau_multi+0.5));
end

figure
plot(tau_multi,gamma_anal_multi)
xlabel('tau_multi [-]')
ylabel('gamma_anal_multi [-]')
figure
subplot(2,1,1)
plot(tau_multi,theta_multi-gamma_num_multi.')
xlabel('tau_multi [-]')
ylabel('theta_multi - gamma_num_multi [-]')
subplot(2,1,2)
plot(tau_multi,theta_multi-gamma_anal_multi)
xlabel('tau_multi [-]')
ylabel('theta_multi - gamma_anal_multi [-]')

theta_multi2=size(tau);
theta_multi2(1:16001)=theta_multi(20000:36000);
theta_multi2(16001:22001)=theta_multi(1:6001);
tau_multi2=size(tau);
tau_multi2(1:16001)=tau_multi(20000:36000);
tau_multi2(16001:22001)=tau_multi(1:6001);
gamma_num_multi2=size(tau);
gamma_num_multi2(1:16001)=gamma_num_multi(20000:36000);
gamma_num_multi2(16001:22001)=gamma_num_multi(1:6001);
gamma_anal_multi2=size(tau);
gamma_anal_multi2(1:16001)=gamma_anal_multi(20000:36000);
gamma_anal_multi2(16001:22001)=gamma_anal_multi(1:6001);

figure 
subplot(2,2,1)
plot(tau,theta-theta_multi2)
xlabel('tau [-]')
ylabel('theta - theta_multi [-]')
subplot(2,2,3)
plot(tau,gamma_num-gamma_num_multi2.')
xlabel('tau [-]')
ylabel('gamma - gamma_num_multi [-]')
subplot(2,2,4)
plot(tau,gamma_num-gamma_anal_multi2.')
xlabel('tau [-]')
ylabel('gamma - gamma_anal_multi [-]')

figure 
plot(tau_multi,gamma_num_multi.' - gamma_anal_multi)
xlabel('tau_multi [-]')
ylabel('gamma_num_multi - gamma_anal_multi [-]')

disp("start Force analysis")
%%%% FORCE ANALYSIS %%%%

force = kf*0.03*(-theta_multi + gamma_num_multi.')./(cos(Pressure_Angle));
figure
plot(tau_multi, force, 'LineWidth', 2)
xlabel('tau_multi [-]')
ylabel('Vibration force')
axis([0, 1, -inf, inf])

figure
plot(theta1, normalForce + force)
xlabel('theta [-]')
ylabel('Normal force + vibration force between follower and cam')
min(normalForce + force);



