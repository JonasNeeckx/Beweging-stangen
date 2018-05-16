function [I_flywheel,meanIndex] = flywheel(instantanious_torque, average_torque,w_mean)

sommation = 0;
sommation_cache = zeros(size(instantanious_torque));

diff_torque = instantanious_torque-average_torque*ones(size(instantanious_torque));
for i = 1:size(instantanious_torque,2)
    sommation = sommation+diff_torque(i);
    sommation_cache(i) = sommation;   
end

[~,minIndex] = min(sommation_cache);
[~,maxIndex] = max(sommation_cache);
[~,meanIndex] = min(abs(sommation_cache));

Amax = sum(diff_torque(minIndex:maxIndex))*2*pi/size(instantanious_torque,2);

I_flywheel = 2*Amax/(1.05*w_mean^2-.95*w_mean^2);