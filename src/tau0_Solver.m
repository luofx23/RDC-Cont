function [tau] = tau0_Solver(uref,u0,beta0,direction)
global nx gamma M
Fx0 = sGjac(u0,beta0,uref);
rho = u0(1:nx); rho_u = u0(nx+1:2*nx); E = u0(2*nx+1:3*nx);rho_z = u0(3*nx+1:4*nx);
p = (gamma-1)*(E-0.5*rho_u.^2./rho);
dfdbeta = - rho_z.*HP(p);
f = [zeros(nx,1);zeros(nx,1);zeros(nx,1);dfdbeta];
F2 = [[M 0*M 0*M 0*M];[0*M M 0*M 0*M];[0*M 0*M M 0*M];[0*M 0*M 0*M M]]*f;
F_lambda0 = [-F2;0];
Ir = speye(4*nx+1,4*nx+1);
z = (Fx0+1e-6*Ir)\(-F_lambda0);
tau_plus = 1/sqrt(z'*z + 1)*[z;1];
tau_neg = -1/sqrt(z'*z + 1)*[z;1];

if direction == 'backward'
    tau = tau_neg;
else
    tau = tau_plus;
end

end

