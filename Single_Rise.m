function Single_Rise(zeta,ks,theta1,Omega_rad,Mass,S,normalForce)

disp(["start Single rise analysis"])
%Rise from 200 till 280 degrees, dwell from 280 till 60 degrees
rise = 200;
dwell = 280;
d_end = 60;
t1 = ((dwell - rise)*pi)/(180*Omega_rad);
tau_end = (d_end-rise+360)/(dwell-rise);
lambda = 0.75/zeta;

kf = Mass*((2*pi*lambda)/(t1))^2 - ks;

tau = 0:(tau_end/22000):tau_end;

numerator = (2*pi*lambda)^2;
denominator = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(numerator, denominator);

approx_theta = (((2*pi)^2*((tau - 1).^3))/factorial(3)) + 1;

Rise = zeros(size(tau));
Rise(1:16001) = S(20000:36000);
Rise(16002:22001) = S(1:6000);
theta = (-30+Rise)./-30;

gamma_num = lsim(sys, theta, tau); %if the system starts from lift and speed equal to zero

% The system doesn't start from lift and speed equal to zero
theta0 = (2*pi)^2 * ((2*zeta - 2*zeta*(4*(zeta^2)-1))/((2*pi*lambda)^3));
theta_dot0 = (2*pi)^2 * ((4*(zeta^2)-1)/((2*pi*lambda)^2));
[A,B,C,D] = tf2ss(numerator,denominator);
X0 = [1/C(2)*theta_dot0; 1/C(2)*theta0];

% compute free response approximation
gamma_approx = lsim(A,B,C,D, theta, tau, X0);

% compare exponential envelopes of numerical and approximate solutions
x_0 = gamma_num(8000) - 1;
lambda_d = lambda * sqrt(1 - zeta^2);
dx_0 = (gamma_num(8001) - gamma_num(7999))/(2*0.0001);
A_num = sqrt(((x_0*2*pi*lambda_d)^2 + (dx_0 + zeta*2*pi*lambda*x_0)^2)/((2*pi*lambda_d)^2));
A_approx = (((2*pi)^2)/((2*pi*lambda)^3))*sqrt(1/(1-zeta^2));

%%%% PLOT SINGLE RISE RESULTS %%%%
interval = 1:22001;     % plotting interval (entire interval of rise and dwell: 1:22001)

% numerical result and input-output difference
figure
subplot(2, 1, 1)
plot(tau(interval), gamma_num(interval), 'LineWidth', 2)
xlabel('tau [-]')
ylabel('numerical solution [-]')
axis([0, tau_end, -inf, inf])
% hold on
% plot(tau(interval), 1 + A_num*exp(-zdeta*2*pi*lambda*(tau(interval) - 1)))

subplot(2, 1, 2)
plot(tau(interval), theta(interval) - gamma_num(interval).', 'LineWidth', 2)
xlabel('tau [-]')
ylabel('input-output difference for numerical solution [-]')
axis([0, tau_end, -inf, inf])

% approximation of response after rise and exponential envelope
figure
subplot(2, 1, 1)
plot(tau(8000:22001), gamma_approx(8000:22001), 'LineWidth', 2)
xlabel('tau [-]')
ylabel('approximation of free response during dwell [-]')
axis([1, tau_end, -inf, inf])
hold on
plot(tau(8000:22001), 1+A_approx*exp(-zeta*2*pi*lambda*(tau(8000:22001) -1)))
hold on
plot(tau(8000:22001),1- A_approx*exp(-zeta*2*pi*lambda*(tau(8000:22001)-1)))
legend('response', 'exponential envelope')

subplot(2, 1, 2)
plot(tau(8000:22001), gamma_approx(8000:22001) - gamma_num(8000:22001), 'LineWidth', 2)
xlabel('tau [-]')
ylabel('difference between numerical and approximate solution [-]')
axis([1, tau_end, -inf, inf])

%%% MULTI RISE ANALYSIS %%%%

% define system
wn = sqrt((ks + kf)/Mass);
tn = (2*pi)/wn;
lambda_tilde = 2/tn;
numerator2 = (2*pi*lambda_tilde)^2;
denominator2 = [1, 2*zeta*(2*pi*lambda_tilde), (2*pi*lambda_tilde)^2];
sys2 = tf(numerator2, denominator2);

% construct full input
tau_MR = theta1/(2*pi);
input_MR1 = 0.015*((tau_MR - 60/360)/(60/360) - sin(2*pi*(tau_MR - 60/360)/(60/360))/(2*pi));
input_MR2 = 0.015 + 0.015*((tau_MR - 120/360)/(60/360) - sin(2*pi*(tau_MR - 120/360)/(60/360))/(2*pi));
input_MR3 = 0.03*(1 - (tau_MR - 200/360)/(80/360) + sin(2*pi*(tau_MR - 200/360)/(80/360))/(2*pi));

input_MR = zeros(size(tau_MR));
input_MR(6001:12001) = input_MR1(6001:12001);
input_MR(12001:18001) = input_MR2(12001:18001);
input_MR(20001:28001) = input_MR3(20001:28001);

% nondimensionalise by dividing by highest peak
input_MR = input_MR/0.03;

% compute fourier series of multi rise input
N = 100; % number of terms
[a, b] = Fseries(tau_MR, input_MR, N);

% compute analytical result with fourier coefficients
c = zeros(N, 1);
d = zeros(N, 1);

for k = 1:N
    c(k) = (-2*zeta*k*b(k)/lambda_tilde + a(k+1)*(1 - ((k/lambda_tilde)^2)))/(((2*zeta*k/lambda_tilde)^2) + ((1 - ((k/lambda_tilde)^2))^2));
    d(k) = (2*zeta*k*a(k+1)/lambda_tilde + b(k)*(1 - ((k/lambda_tilde)^2)))/(((2*zeta*k/lambda_tilde)^2) + ((1 - ((k/lambda_tilde)^2))^2));
end

gamma_anal_MR = 0.5*a(1)*ones(size(tau_MR));

for k = 1:N
    gamma_anal_MR = gamma_anal_MR + c(k)*cos(2*pi*k*(tau_MR)) + d(k)*sin(2*pi*k*(tau_MR));
end

% compute multi rise numerical response with lsim
gamma_num_MR = lsim(sys2, input_MR, tau_MR);

%%%% PLOT MULTI RISE RESULTS %%%%
figure
subplot(2, 1, 1)
plot(tau_MR, gamma_num_MR, 'LineWidth', 2)
xlabel('tau [-]')
ylabel('numerical solution [-]')
axis([0, 1, -inf, inf])

subplot(2, 1, 2)
plot(tau_MR, input_MR.' - gamma_num_MR, 'LineWidth', 2)
xlabel('tau [-]')
ylabel('input-output difference for numerical solution [-]')
axis([0, 1, -inf, inf])

figure
plot(tau_MR, gamma_num_MR - gamma_anal_MR.', 'LineWidth', 2)
xlabel('tau [-]')
ylabel(strcat('difference between numerical and analytical solution with ', num2str(N), ' fourier terms [-]'))
axis([0, 1, -inf, inf])

%%%% COMPARE SINGLE AND MULTI RISE %%%%
figure
subplot(2, 1, 1)
plot(tau_MR(20001:end), -theta(1:16000) - input_MR(20001:end) + 1, 'LineWidth', 2)
xlabel('tau [-]')
ylabel('difference in input between multi-rise and single-rise [-]')
axis([120/360, 1, -inf, inf])

subplot(2, 1, 2)
plot(tau_MR(20001:end), -gamma_num(1:16000) - gamma_num_MR(20001:end) + 1, 'LineWidth', 2)
xlabel('tau [-]')
ylabel('difference in response between multi-rise and single-rise [-]')
axis([120/360, 1, -inf, inf])

%%%% FORCE ANALYSIS %%%%

force = kf*0.03*(input_MR - gamma_num_MR.');
figure
plot(tau_MR, force, 'LineWidth', 2)
xlabel('tau [-]')
ylabel('resulting force between follower and cam')
axis([0, 1, -inf, inf])

figure
plot(theta1, normalForce + force)



