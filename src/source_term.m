function [f] = source_term(u,beta,Da,Tc)
global gamma
rho = u(1:end/4); rho_u = u(end/4+1:end/2); E = u(end/2+1:3/4*end);rho_z = u(3/4*end+1:end); 
alpha = sqrt(gamma)*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
Ap = 0.2; Am = 1.0; q = 24.6; Tvn = 5.8; Ea = 11.5;
p = (gamma-1)*(E-0.5*rho_u.^2./rho);
T = p./rho;
omega = (T>Tc).*Da.*(rho - rho_z).*exp(-Ea*(1./T - 1./Tvn));
f1 = alpha*(Ap*HP(p)-Am*sqrt(rho.*p));
f2 = zeros(length(f1),1);
f3 = alpha/(gamma-1)*(Ap*HP(p)-Am*T.*sqrt(rho.*p)) + omega*q;
f4 = omega - rho_z.*beta.*HP(p) + alpha*(Ap*HP(p)-Am*sqrt(rho.*p)).*rho_z./rho;
f = [f1;f2;f3;f4];
end

