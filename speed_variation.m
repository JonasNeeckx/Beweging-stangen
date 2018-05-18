function w = speed_variation(average_torque,instantaneous_torque,I, i_mean, w_mean) 
%SPEED_VARIATION calculates the variation in rotational speed of the cam.
%This is supposed to be smaller than 5% compared to the mean rotational
%speed. 
%
%INPUT
%@param average_torque
%           The average torque posed upon the cam. A scalar is expected [Nm]
%@param Instantaneous_torque
%           The torque applied on the cam at a discrete moment. A vector is
%           expected. [Nm]
%@param I
%           The moment of inertia of the cam and flywheel together. [kgm^2]
%@param i_mean
%           The index where the average rotational speed is reached
%@param w_mean
%           The average rotational speed [rad/s]
%
%OUTPUT
%@result w
%           The instantaneous rotational speed [rad/s]
disp(["Analysis of the speed variation"])

theta = 0:0.01:359.99;
dim = size(instantaneous_torque,2);
w = zeros(dim,1);
w(i_mean) = w_mean;
for i = i_mean+1:i_mean+dim-1
    if i > dim
        index = i-dim;
    else
        index=i;
    end
    
    if index == 1
        prev = dim;
    else
        prev = index-1;
    end
    
    w(index) = w(prev)-2*pi/size(instantaneous_torque,2)*(instantaneous_torque(index)-average_torque)/(I*w(prev));
end

figure 
subplot(211)
plot(theta, w)
xlabel('theta [°]')
ylabel('Rotational speed [rad/s]')

subplot(212)
plot(theta, w/w_mean)
xlabel('theta [°]')
ylabel('Relative variation of the rotational speed[~]')