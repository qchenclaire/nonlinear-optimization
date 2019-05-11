function [p, iters, flag]=steihaug_CG(B, g,radius, tol)
    n=length(g);
    p=zeros(n,1);
    iters = 0;
    r=g;
    s=-g;
    term = tol * norm(g);
    while norm(r) > term
        iters = iters+1;
        % solve tau
        a=s'*s;
        b=2*p'*s;
        c=p'*p-radius*radius;
        para = [a b c];
        soltau = roots(para);
        if soltau(1)>0
            tau = soltau(1);
        else
            tau = soltau(2);
        end
        tmp =  s'*B*s;
        if tmp > 0
            alpha = r'*r/tmp;
        else        
            p=p+tau*s;
            flag = -1;
            return
        end
        if norm(p+alpha*s)<radius
            p=p+alpha*s;
        else
            p=p+tau*s;
            flag = 1;
            return;
        end
        denom = r'*r;
        r=r+alpha*B*s;
        beta = r'*r/denom;
        s=-r+beta*s;
    end
    flag = 0;
    return;
end

             
            
            