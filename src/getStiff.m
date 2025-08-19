function [K,Kx,M] = getStiff(L,nx)
dx = L/nx;
v1 = -0.5*ones(nx-1,1);
Kx = spdiags(v1,-1,nx,nx) + (spdiags(-v1,-1,nx,nx)');
Kx(1,end) = -0.5;
Kx(end,1) = 0.5;

v2 = 2*ones(nx,1);
v3 = -1*ones(nx-1,1);
K = spdiags(v2,0,nx,nx) + spdiags(v3,-1,nx,nx) + (spdiags(v3,-1,nx,nx)');
K(1,end) = -1;
K(end,1) = -1;
K = K/dx;

v4 = 4*ones(nx,1);
v5 = 1*ones(nx-1,1);
M = spdiags(v4,0,nx,nx) + spdiags(v5,-1,nx,nx) + (spdiags(v5,-1,nx,nx)');
M(1,end) = 1;
M(end,1) = 1;
M = M/6*dx;
end

