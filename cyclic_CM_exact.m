function [x, F, J,iter, status]=cyclic_CM_exact(x0, H, maxit)
n = length(x0);
th = max([norm(H*x0),1])* 10^(-6);
x = x0;
status = 1;
iter = 0;
for k=1:maxit 
    J = H*x;
    if norm(J) <= th 
        status = 0;
        break; 
    end
    iter = iter + 1;
    i = mod(k, n)+1;
    Hi = H(i,:);
    Hii = Hi(i);
    alpha = Hi*x/Hii;
    x(i) = x(i) - alpha;
    
    F = 0.5*x'*H*x;
    g=sprintf('%f ', x);
    fprintf('iters:%i x: %s F:%s Norm(J):%f\n', k, g,F, norm(J));  
end
end
