function [rho,u,p,z,s] = state_reconstruct(sol)
global gamma
nx = round((length(sol)-1)/4);
rho = sol(1:nx);
u = sol(nx+1:2*nx)./rho;
p = (sol(2*nx+1:3*nx)-0.5.*rho.*u.^2)*(gamma-1);
z = sol(3*nx+1:4*nx)./rho;
s = sol(end);
end

