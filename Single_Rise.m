function Single_Rise(zeta,Springconstant_optimal,theta,Omega_rad,Mass,S)

disp(["start Single rise analysis"])
%Rise from 200 till 280 degrees, dwell from 280 till 60 degrees
rise = 200;
dwell = 280;
d_end = 60;
t1 = ((dwell - rise)*pi)/(180*Omega_rad);
tau_end = (d_end-rise+360)/(dwell-rise);
lambda = 0.75/zeta;

kfmin = Mass*((0.75*2*pi)/(zeta*t1))^2 - Springconstant_optimal;

wn = sqrt((Springconstant_optimal + kfmin)/Mass);
tn = (2*pi)/wn;
lambda2 = t1/tn;

decay1 = exp(-zeta*2*pi*lambda*(tau_end - 1));
decay2 = exp(-zeta*2*pi*lambda);


% find indices of events
rise_index = 100*rise;
dwell_index = 100*dwell;
d_end_index = 36000 + 100*d_end;

% create tau for single rise by extending theta vector (dwell ends after
% vector ends)
extended_time = [theta, (theta(2:36000) + theta(36000))];
rescaled_time = (extended_time - extended_time(rise_index))/(extended_time(dwell_index) - extended_time(rise_index));
tau = rescaled_time(rise_index:d_end_index);

numerator = (2*pi*lambda)^2;
denominator = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(numerator, denominator);

approx_input = (((2*pi)^2*((tau - 1).^3))/factorial(3)) + 1;

extended_rise = [S, S(2:36000)];
rescaled_rise = (extended_rise - extended_rise(rise_index))/(extended_rise(dwell_index) - extended_rise(rise_index));
theta = rescaled_rise(rise_index:d_end_index);

gamma_num = lsim(sys, theta, tau);

% starting conditions for free response approximation
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
% plot(tau_SR(interval), 1 + A_num*exp(-zdeta*2*pi*lambda*(tau_SR(interval) - 1)))

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

