function [out] = HP(p)
global gamma
Y = (1+(gamma-1)/2)^(-gamma/(gamma-1));
tm1 = heaviside(1 - p);
tm2 = (1 - heaviside(p - Y).*((p - Y)./(1 - Y)));
out = tm1.*tm2;
end

