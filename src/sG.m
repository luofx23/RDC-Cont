function r = sG(u,beta,uref)

nx = round((size(u,1)-1)/4);nu = 0.0075;
global Ks Kx M gamma
s = u(end);
rho = u(1:nx); rho_u = u(nx+1:2*nx); E = u(2*nx+1:3*nx);rho_z = u(3*nx+1:4*nx);
p = (gamma-1)*(E-0.5*rho_u.^2./rho);
K=[[nu*Ks-s*Kx 0*Kx 0*Kx 0*Kx];
   [0*Kx nu*Ks-s*Kx 0*Kx 0*Kx];
   [0*Kx 0*Kx nu*Ks-s*Kx 0*Kx];
   [0*Kx 0*Kx 0*Kx nu*Ks-s*Kx]]; 
F1 = [-Kx*(rho_u); -Kx*(rho_u.^2./rho+p); -Kx*(rho_u./rho.*(E+p)); -Kx*(rho_u./rho.*rho_z)]; 
f = nodalf(u(1:4*nx),beta);
F2 = [[M 0*M 0*M 0*M];[0*M M 0*M 0*M];[0*M 0*M M 0*M];[0*M 0*M 0*M M]]*f;
phase_condition_u = u(1:nx)'*(Kx*uref(1:nx))/norm(Kx*uref(1:nx));
r = [K*u(1:4*nx) - (F1 + F2);phase_condition_u]; 