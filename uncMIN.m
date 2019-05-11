function [x, F, G, H, iter, status] = uncMIN(fun, x0, step, maxit, printlevel, tol)
%Input
% fun -- type string  that holds the name of a Matlab m-function
% x0 -- an intial guess at a zero
% step -- if 0, steepest-descent; modified-Newton otherwise
% maxit -- the maximum number of iterations allowed
% printlevel -- the amount of printout required
% tol -- final stopping tolerance

%Output
%x, F, G, H -- final iterate, function value, Gradient vector, Hessian
%matrix
%iter -- total number of iterations performed
%status -- 0 if final stopping tolerance was obtained; 1 otherwise
x = x0;
status = 1;
iter = 0;
[F0, G0, H0] = feval(fun, x0);
for i=1:maxit   
    [F, G, H] = feval(fun, x);
    if printlevel ~= 0
        g=sprintf('%f ', x);
        k=sprintf('%f ', F);
        fprintf('iter:%i x: %s F:%s Norm(F):%f\n', i,g,k, norm(F));
    end
    iter = iter + 1;
    if norm(F)  <  max(norm(F0),1) * tol
        status = 0;
        break;
    end
    if step ==0
        p = -G;
    else
        [B, flag] = modNewton(H, 10^7);
        p = B\(-G);
    end
    % backtracking -Armijo linesearch
    alpha = 1;
    eta = 0.5;
    tao = 0.5;
    while 1
        [F2, G2, H2] = feval(fun, x + alpha * p);
        if (F2 > F + eta * alpha * G' * p)
            alpha = alpha * tao;
        else
            break;
        end
    end
    x = x + alpha * p;    
end
    if norm(F)  <  max(norm(F0),1) * tol        
        status = 0;
    end   
end