function [x, F, J,iter, status]=GS_CM_fixed(x0, H, maxit)
n = length(x0);
alpha = 1/norm(H,'fro');
th = max(norm(H*x0),1)* 10^(-6);
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
    [~, i] = max(abs(J));
    x(i) = x(i) - alpha * J(i);   
    F = 0.5*x'*H*x;
    g=sprintf('%f ', x);
    fprintf('iters:%i x: %s F:%s Norm(J):%f\n', k, g,F, norm(J));    
end
end