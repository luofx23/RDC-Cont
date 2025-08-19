function [out] = HP_dot(p)
global gamma
Y = (1+(gamma-1)/2)^(-gamma/(gamma-1));
out = heaviside(1-p).*(-(1/(1-Y))*heaviside(p-Y));
end

