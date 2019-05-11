
figure;
for y=0:0.02:1.5
syms x
eqn = f(x)== (y*y);
solx = vpasolve(eqn,x,[0 Inf]);
if solx
    u = -1/(1+double(solx));
    v = -1/(3+double(solx));
else
    u = -1;
    v = -1/3;
end
plot(u, v, '.');
hold on
circles(0,0,y,'facealpha',0);
hold on
end