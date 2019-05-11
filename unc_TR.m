function [x, F, G, H,iter, status]=unc_TR(fun, x0, maxit, printlevel,tol)
gamma_d = 0.5;
gamma_i = 1.5;
eta_s = 0.001;
eta_vs = 0.9;
radius=1;
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
    
    [s, iters, flag]=steihaug_CG(H, G,radius, tol);
    [F_s,~, ~] = feval(fun, x+s);
    thou = (F-F_s)/(-G'*s-0.5*s'*H*s);
    if thou >= eta_vs
        x = x+s;
        radius = radius * gamma_i;
    else
        if thou >= eta_s
            x = x+s;
        else
            radius = radius * gamma_d;
        end
    end
    
end
if norm(F)  <  max(norm(F0),1) * tol        
    status = 0;
end   
end