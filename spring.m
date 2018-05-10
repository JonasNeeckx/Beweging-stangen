function [Springconstant_optimal, Fv0_optimal] = spring(Lift, F_load, F_inert, Pressure_Angle)
optimal_N = inf;

for Fv0 = 50:0.01:300
    Springconstant= (331.7 +Fv0 + 0.008841)/14.99;
    N = (F_load + Fv0 + Springconstant*Lift + F_inert)/cos(Pressure_Angle);
    Nmax=max(N);
    Nmin=min(N);
    if (Nmax < optimal_N) && (Nmin >= 0)
        optimal_N = N;
        Springconstant_optimal = Springconstant;
        Fv0_optimal = Fv0;
    end
end