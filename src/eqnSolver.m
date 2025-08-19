function [solution,tol] = eqnSolver(initial_guess,beta)
x = initial_guess;
nx = round((size(x,1)-1)/4);
Ir = speye(4*nx+1,4*nx+1);
for i = 1:500
    Jac = sGjac(x,beta,initial_guess);
    dJ = Jac;
    r = sG(x,beta,initial_guess);
    x_next = x - (dJ + 1e-5*Ir)\r;
    x = x_next;
    if norm(sG(x_next,beta,initial_guess),'inf')<1e-6
        break
    end
end
solution = x;
tol = norm(sG(x_next,beta,initial_guess),'inf');
end

