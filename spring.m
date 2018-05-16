function [Springconstant_optimal, Fv0_optimal,optimal_N] = spring(S, F_load, F_inert, Pressure_Angle)
optimal_N = inf;

for Fv0 = 200:.1:300
    for Springconstant = 1.5:0.01:3
        N = F_load +  + F_inert + (Fv0*ones(size(F_load)) + Springconstant*S)./cos(Pressure_Angle);
        Nmax=max(N);
        Nmin=min(N);
        if (Nmax < optimal_N) && (Nmin >= 0)
            optimal_N = Nmax;
            Springconstant_optimal = Springconstant;
            Fv0_optimal = Fv0;
        end
    end
end
