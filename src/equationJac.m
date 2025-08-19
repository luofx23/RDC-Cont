function [Gu] = equationJac(u,s,beta)
nu = 0.0075;
global Ks Kx M nx gamma dx
rho = u(1:nx); rho_u = u(nx+1:2*nx); E = u(2*nx+1:3*nx);rho_z = u(3*nx+1:4*nx);
[f11,f12,f13,f14,f21,f22,f23,f24,f31,f32,f33,f34,f41,f42,f43,f44] = njac(u(1:4*nx),beta);

K=[[nu*Ks-s*Kx 0*Kx 0*Kx 0*Kx];
   [0*Kx nu*Ks-s*Kx 0*Kx 0*Kx];
   [0*Kx 0*Kx nu*Ks-s*Kx 0*Kx];
   [0*Kx 0*Kx 0*Kx nu*Ks-s*Kx]];
F12 = -Kx*spdiags(ones(nx,1),0,nx,nx);

F21 = -Kx*spdiags((0.5*(gamma-1)-1)*rho_u.^2./(rho.^2),0,nx,nx);
F22 = -Kx*spdiags(2*(1-0.5*(gamma-1))*rho_u./rho,0,nx,nx);
F23 = -Kx*spdiags((gamma-1)*ones(nx,1),0,nx,nx);
F24 = 0*F23;

F31 = -Kx*spdiags(0.5*(gamma-1)*rho_u.^2./(rho.^2)-gamma*E.*rho_u./(rho.^2),0,nx,nx);
F32 = -Kx*spdiags(gamma*E./rho - (gamma-1)*rho_u./rho,0,nx,nx);
F33 = -Kx*spdiags(gamma*rho_u./rho,0,nx,nx);
F34 = 0*F33;

F41 = -Kx*spdiags(-rho_z.*rho_u./(rho.^2),0,nx,nx);
F42 = -Kx*spdiags(rho_z./rho,0,nx,nx);
F43 = 0*F42;
F44 = -Kx*spdiags(rho_u./rho,0,nx,nx);
 
F1u=[[0*F12 F12 0*F12 0*F12];
     [F21 F22 F23 F24];
     [F31 F32 F33 F34];
     [F41 F42 F43 F44]];
F2u=[[M 0*M 0*M 0*M];[0*M M 0*M 0*M];[0*M 0*M M 0*M];[0*M 0*M 0*M M]]*...
    [[spdiags(f11,0,nx,nx) spdiags(f12,0,nx,nx) spdiags(f13,0,nx,nx) spdiags(f14,0,nx,nx)];
     [0*spdiags(f21,0,nx,nx) 0*spdiags(f22,0,nx,nx) 0*spdiags(f23,0,nx,nx) 0*spdiags(f24,0,nx,nx)];
     [spdiags(f31,0,nx,nx) spdiags(f32,0,nx,nx) spdiags(f33,0,nx,nx) spdiags(f34,0,nx,nx)];
     [spdiags(f41,0,nx,nx) spdiags(f42,0,nx,nx) spdiags(f43,0,nx,nx) spdiags(f44,0,nx,nx)]];
Gu = K-(F1u+F2u);
Gu = -Gu/dx;

%[row,col,v] = find(Gu);
%row = row(abs(v)>1e-10);
%col = col(abs(v)>1e-10);
%v = v(abs(v)>1e-10);
%Gu = sparse(row,col,v,9600,9600);
end

