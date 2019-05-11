function [x, F, J, iter, status] = Newton( Fun, x0, maxit, printlevel, tol)
%Input
% Fun -- type string  that holds the name of a Matlab m-function
% x0 -- an intial guess at a zero
% maxit -- the maximum number of iterations allowed
% printlevel -- the amount of printout required
% tol -- final stopping tolerance

%Output
%x, F, J -- final iterate, function value, Jacobian matrix
%iter -- total number of iterations performed
%status -- 0 if final stopping tolerance was obtained; 1 otherwise
x = x0;
status = 1;
[F0, J] = feval(Fun, x0);
for i=1:maxit   
    [F, J] = feval(Fun, x);
    if printlevel ~= 0
        g=sprintf('%f ', x);
        k=sprintf('%f ', F);
        fprintf('iter:%i x: %s F:%s Norm(F):%f\n', i,g,k, norm(F));
    end
    if norm(F)  <  norm(F0) * tol
        iter = i-1;
        status = 0;
        break;
    end    
    s = J\(-F);
    x = x + s;    
end
    if norm(F)  <  norm(F0) * tol        
        status = 0;
    end
    iter = maxit;
end

    
    