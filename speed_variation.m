function w = speed_variation(average_torque,instantanious_torque,I, wstart) 

w = zeros(size(instantanious_torque));
w(1) = wstart;
for i = 2:size(instantanious_torque,2)
    w(i) = w(i-1)+2*pi/size(instantanious_torque,2)*(instantanious_torque(i)-average_torque)/(I*w(i-1));
end