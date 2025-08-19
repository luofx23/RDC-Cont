function [f11,f12,f13,f14,f21,f22,f23,f24,f31,f32,f33,f34,f41,f42,f43,f44]=njac(u,beta) % local (no spat.derivatives) Jacobian 
global gamma
rho = u(1:end/4); rho_u = u(end/4+1:end/2); E = u(end/2+1:3/4*end);rho_z = u(3/4*end+1:end);
v = rho_u./rho;z = rho_z./rho;
alpha = sqrt(gamma)*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
Ap = 0.2; Am = 1; q = 24.6; Da = 128.7; Tvn = 5.8; Ea = 11.5;
p = (gamma-1)*(E-0.5*rho_u.^2./rho);
T = p./rho;
omega = Da*(rho - rho_z).*exp(-Ea*(1./T - 1./Tvn));
f1 = alpha*(Ap*HP(p)-Am*sqrt(rho.*p));

f11 = alpha*(Ap*HP_dot(p)*0.5.*(gamma-1).*v.^2 - 0.5*Am*(sqrt(rho./p)*0.5.*(gamma-1).*v.^2 + sqrt(p./rho)));
f12 = alpha*(Ap*HP_dot(p)*(1-gamma).*v - 0.5*Am*sqrt(rho./p)*(1-gamma).*v);
f13 = alpha*(Ap*HP_dot(p)*(gamma-1) - 0.5*Am*sqrt(rho./p)*(gamma-1));
f14 = zeros(length(f11),1);

f21 = zeros(length(f11),1);f22 = zeros(length(f11),1);
f23 = zeros(length(f11),1);f24 = zeros(length(f11),1);

dTdrho = 1./(rho*1)*(gamma-1)*0.5.*v.^2 - p./(rho.^2);
dTdrho_u = 1./rho.*(1-gamma).*v;
dTdE = 1./rho.*(gamma-1);

dwdrho = Da*exp(-Ea*(1./T - 1./Tvn)) + 1./(T.^2)*Ea.*omega.*dTdrho;
dwdrho_u = 1./(T.^2)*Ea.*omega.*dTdrho_u;
dwdE = 1./(T.^2)*Ea.*omega.*dTdE;
dwdrho_z = -Da*exp(-Ea*(1./T - 1./Tvn));

f31 = alpha/(gamma-1).*(Ap*HP_dot(p)*(gamma-1)*0.5.*v.^2 - 0.5*Am*T.*(sqrt(rho./p)*(gamma-1)*0.5.*v.^2 + sqrt(p./rho))...
                        -Am*sqrt(p.*rho).*dTdrho) + dwdrho*q;
f32 = alpha/(gamma-1).*(Ap*HP_dot(p)*(1-gamma).*v - 0.5*Am*T.*sqrt(rho./p)*(1-gamma).*v...
                        -Am*sqrt(p.*rho).*dTdrho_u) + dwdrho_u*q;
f33 = alpha/(gamma-1).*(Ap*HP_dot(p)*(gamma-1) - 0.5*Am*T.*sqrt(rho./p)*(gamma-1)...
                        -Am*sqrt(p.*rho).*dTdE) + dwdE*q;
f34 = dwdrho_z*q;

f41 = dwdrho - beta*rho_z.*HP_dot(p)*(gamma-1)*0.5.*v.^2 + f11.*z - z./rho.*f1;
f42 = dwdrho_u - beta*rho_z.*HP_dot(p)*(1-gamma).*v + f12.*z;
f43 = dwdE - beta*rho_z.*HP_dot(p)*(gamma-1) + f13.*z;
f44 = dwdrho_z - beta*HP(p) + 1./rho.*f1;

end
