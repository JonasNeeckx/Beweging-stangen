function [Springconstant_optimal, Fv0_optimal,optimal_Nmax] = spring(S, F_load, F_inert, Pressure_Angle)
%SPRING calculates the optimal values for the spring applied on a cam.
%These are obtained by looking for an arrangement where the maximal normal
%force is minimal. 
%
%INPUT
%@param S
%           The lift of the cam. A column is expected. [mm]
%@param F_load
%           The external load applied on the cam normal to the cam's
%           surface. A column of different values is expected. [N]
%@param F_intert
%           The inertial force of the follower normal to the surface. A column is expected [N]
%@param Pressure_Angle
%           The angle between the normal on the cam's surface and the line
%           connecting the follower with the rotation point. [rad]
%
%OUTPUT
%@result Springconstant_optimal
%           The optimal spring constant [N/mm]
%@result Fv0_optimal
%           The optimal prestress [N]
disp("calculating optimal spring setting")

optimal_Nmax = inf;
optimal_Nmean = inf; 

for Fv0 = 220:0.1:230
    for Springconstant = 5:0.01:5.5
        N = F_load +  4*F_inert + (Fv0*ones(size(F_load)) + Springconstant*S)./cos(Pressure_Angle);
        Nmax=max(N);
        Nmin=min(N);
        Nmean = mean(N);
        if ((Nmax < optimal_Nmax) && (Nmin >= 0)||((Nmax == optimal_Nmax) && (Nmean < optimal_Nmean)))
            optimal_Nmax = Nmax;
            optimal_Nmean = Nmean;
            Springconstant_optimal = Springconstant;
            Fv0_optimal = Fv0;
        end
    end
end
