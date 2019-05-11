function [ B, flag] = modNewton(H, beta)
% input
% H symmetric
% beta > 1, upper bound on the required condition number of modified matrix

% output
% B --positive-definite matrix 
% flag --0 if no modification; 1 otherwise

if H==0
    epsilon = 1;
else
    epsilon = norm(H)/beta;
end
[V, D] = eig(H);
D_bar = D + diag(max(max(0,-2*diag(D)),epsilon - diag(D)));
B = V*D_bar*V';
flag = ~isequal(D, D_bar);
