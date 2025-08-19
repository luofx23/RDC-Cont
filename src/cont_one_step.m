function [u,beta,tau,tol] = cont_one_step(tau0,u0,beta0,uref,ds)
%CONT_ONE_STEP 此处显示有关此函数的摘要
%   此处显示详细说明
Fx0 = sGjac(u0,beta0,uref);
global nx M gamma
%s = u0(end);
u_whole = u0;%impose phase constrain on state(u(0) = 1)
rho = u_whole(1:nx); rho_u = u_whole(nx+1:2*nx); E = u_whole(2*nx+1:3*nx);rho_z = u_whole(3*nx+1:4*nx);
p = (gamma-1)*(E-0.5*rho_u.^2./rho);
%omega = Da*(rho - rho_z).*exp(-Ea*(1./T - 1./Tvn));
dfdbeta = - rho_z.*HP(p);
f = [zeros(nx,1);zeros(nx,1);zeros(nx,1);dfdbeta];
F2 = [[M 0*M 0*M 0*M];[0*M M 0*M 0*M];[0*M 0*M M 0*M];[0*M 0*M 0*M M]]*f;
F_lambda0 = [-F2;0];
Ir = speye(4*nx+1,4*nx+1);
z = (Fx0+1e-6*Ir)\(-F_lambda0);
tau_plus = 1/sqrt(z'*z + 1)*[z;1];
tau_neg = -1/sqrt(z'*z + 1)*[z;1];
if tau_plus'*tau0>0
    tau = tau_plus;
else
    tau = tau_neg;
end
%[~,~,tau] = svds([Fx0 F_lambda0],1,'smallest');

%Solve H equations with Newton Iterations
x = [u0;beta0];
mu = 1e-5;
I = speye(4*nx+2,4*nx+2);
for i = 1:500
    [r,Jac] = sH(u0,beta0,tau,uref,ds,x(1:end-1),x(end));
    x_next = x - (Jac + mu*I)\(r);
    x = x_next;
    if norm(sH(u0,beta0,tau,uref,ds,x(1:end-1),x(end)),'inf')<3e-5
        break
    end
    if sum(isnan(real(x)))>0
        disp('Newton iteration Failed. Please change ds');
        break
    end
end
solution = x;
u = solution(1:end-1);beta = solution(end);tol = norm(sH(u0,beta0,tau,uref,ds,x(1:end-1),x(end)),'inf');
end

