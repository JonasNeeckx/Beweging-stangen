function [Springconstant_optimal, Fv0_optimal,optimal_N] = spring(S, F_load, F_inert, Pressure_Angle)
optimal_N = inf;

for Fv0 = 50:0.01:300
    for Springconstant = 0:0.01:40
    N = (F_load + Fv0 + Springconstant*S + F_inert)/cos(Pressure_Angle);
    Nmax=max(N);
    Nmin=min(N);
    if (Nmax < optimal_N) && (Nmin >= 0)
        optimal_N = N;
        Springconstant_optimal = Springconstant;
        Fv0_optimal = Fv0;
    end
    end
end