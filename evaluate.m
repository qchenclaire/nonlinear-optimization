n=1000;
runs = [];
runs_status = [];
a = [];
b = [];
d = [];
a_status=[];
b_status=[];
d_status=[];
c = [];
c_status=[]; 
for cn=[10,100,1000,10000]
       
    rc = 1/cn;
    H = sprandsym(n, 0.2, rc,1);
    x0=-1000 + (1000+1000)*rand(n,1);
    [x, F, J,iter, status]=cyclic_CM_exact(x0, H, 65000000);
    a = [a iter];
    a_status = [a_status status];
    [x, F, J,iter, status]=cyclic_CM_fixed(x0, H, 65000000);
    b = [b iter];
    b_status = [b_status status];
    [x, F, J,iter, status]=GS_CM_fixed(x0, H, 65000000);
    d = [d iter];
    d_status = [d_status status];
    
    
   
    [x, F, J,iter, status]=random_CM_fixed(x0, H, 65000000);
    c = [c iter];
    c_status = [c_status status];        
    
    %runs = [runs ; c];
    %runs_status = [runs_status ; c_status];
end

