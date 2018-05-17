function plot_power(average_power, instantaneous_power, instantaneous_power_no_e,theta)
%CONTROL_POWER plots the power and the control of the power calculation
%
%INPUT
%@param average_power
%           The average power [W]
%@param instantaneous_power
%           The power of a normal cam at a discrete moment [W]
%@param instantaneous_power_no_e
%           The power of a cam without excentricity at a discrete moment
%           [W]
%@param theta
%           The angle where each value is important [°]
disp(['Power calculated'])

control_power = instantaneous_power-instantaneous_power_no_e;

figure 
subplot(211)
plot(theta, instantaneous_power, theta, average_power*ones(size(instantaneous_power)))
xlabel('theta [°]')
ylabel('Power [W]')
legend('Instantaneous power', 'Average power')

subplot(212)
plot(theta, control_power)
xlabel('theta [°]')
ylabel('Difference between the power with and without excentricity [W]')