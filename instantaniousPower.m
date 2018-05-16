function instantaniousPower = instantaniousPower(N_tot,alpha,rpm,V)
%INSTANTANIOUSPOWER calculates the instantanious power needed to drive the
%cam. 
%
%INPUT
%@param N_tot
%           the total force normal on the cam surface
%@param alpha
%           the pressure angle on the cam
%@param rpm
%           the rotational speed in rotations per minute
%@param V
%           the change in elevation of the cam
%
%OUTPUT
%@result instantaniousPower
%           the instantanious power needed to drive the cam. 

omega = 2*pi*rpm/60;
instantaniousPower= zeros(size(V));

%calculate the torque
for i = 1:size(V,2)
    instantaniousPower(i) = N_tot(i)*cos(alpha(i))*V(i)*omega*10^-3;
end


