function [I_flywheel,meanIndex] = flywheel(instantaneous_torque, average_torque,w_mean)
%FLYWHEEL calculates the needed inertia of the flywheel in order to have a
%speed variation beneath 5% for a cam. 
%
%INPUT
%@param instaneous_torque
%           The torque applied upon the cam at a discrete moment. A vector
%           is expected. [Nm]
%@param average_torque
%           The average torque on the cam. A scalar is expected. [Nm]
%@param w_mean
%           The average rotational speed of the cam in rad/s
%
%OUTPUT
%@result I_flywheel
%           The inertia necessary in order to obtain a speed variation that
%           is sufficiently low. 
%@result meanIndex
%           The index where the mean velocity is reached, used in further
%           calculations. 
disp(["Calculating the flyweel settings"])

dim = size(instantaneous_torque,2);
sommation = 0;
sommation_cache = zeros(1,dim);

diff_torque = instantaneous_torque-average_torque*ones(1,dim);
for i = 1:dim
    sommation = sommation+diff_torque(i);
    sommation_cache(i) = sommation;   
end

[~,minIndex] = min(sommation_cache);
[~,maxIndex] = max(sommation_cache);
[~,meanIndex] = min(abs(sommation_cache-mean(sommation_cache)));

Amax = sum(diff_torque(minIndex:maxIndex))*2*pi/dim;
K = 0.06;

I_flywheel = Amax/(K*w_mean^2);