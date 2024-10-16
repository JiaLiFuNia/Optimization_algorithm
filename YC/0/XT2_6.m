clear;
close all;
beta = 0.5;
sigma = 0.2;
xk = [-1,1]';
dk = [1,1]';
m = 0;
M = 20;

syms x y;
fx = 100*(x^2-y)^2+(x-1)^2;
gx = [diff(fx,x),diff(fx,y)];

while m <= M
    if(subs(fx,[x;y],xk+beta^m*dk) <= (subs(fx,[x;y],xk) + sigma*beta^m*subs(gx,[x;y],xk)*dk))
        mk = m;
        break;
    else
        m = m + 1;
    end
end
xk = xk + beta^mk*dk;
vpa(subs(fx,[x;y],xk),6)