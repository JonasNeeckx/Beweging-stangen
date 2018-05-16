function e_optimal = optimalExcentricity(S, V, R0)
%OPTIMALEXCENTRICITY calculates the optimal excentricity for the cam by
%examining the maximal pressure angle at every excentricity. As a check the
%optimal excentricity should also be the excentricity where the maximal and
%minimal pressure angle have the same amplitude. 
%
%INPUT
%@param S
%           S describes the elevation on the cam
%@param V
%           V describes the change in elevation on the cam
%@param R0
%           R0 is the pitch circle radius
%
%OUTPUT
%@result e_optimal
%           e_optimal is the optimal excentricity value

optimal_alpha = inf;
alpha = zeros(size(S));

for e = -10:0.01:-4 
    disp(e);
    for i = 1:size(S, 2)
        alpha(i) = atan((V(i)-e)/(sqrt(R0^2-e^2)+S(i)));
    end
    max_alpha = max(abs(alpha));
    if max_alpha < optimal_alpha
        e_optimal = e;
        optimal_alpha = max_alpha;
        
    end
end