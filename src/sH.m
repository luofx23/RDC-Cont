function [r,Jac] = sH(u0,beta0,tau,uref,ds,u,beta)
global gamma M nx
F = sG(u,beta,uref);
Fx = sGjac(u,beta,uref);
u_whole = u;
rho = u_whole(1:nx); rho_u = u_whole(nx+1:2*nx); E = u_whole(2*nx+1:3*nx);rho_z = u_whole(3*nx+1:4*nx);
p = (gamma-1)*(E-0.5*rho_u.^2./rho);
dfdbeta = - rho_z.*HP(p);
f = [zeros(nx,1);zeros(nx,1);zeros(nx,1);dfdbeta];
F2 = [[M 0*M 0*M 0*M];[0*M M 0*M 0*M];[0*M 0*M M 0*M];[0*M 0*M 0*M M]]*f;
F_lambda = [-F2;0];

%Get tau as tangent vector
tau_x0 = tau(1:end-1);tau_lm0 = tau(end);

r = [F;tau'*[u-u0;beta-beta0]-ds];
Jac = [[Fx F_lambda];[tau_x0' tau_lm0]];
end

