%% 基本牛顿法
clear;
close all;
format long;

eps = 1e-4;

for n = [10,100,1000,2000,5000,10000]

    x0 = 0.01*ones(n,1);
    fx = @(x) fun(x,n);
    gfx = @(x) gradient1(x,n);
    Gfx = @(x) hessian1(x,n);

    k = 0;
    nd_x = x0;
    while(1)
        gk = gfx(nd_x);
        if(norm(gk) <= eps)
            break;
        end
        dk = lsqminnorm(Gfx(nd_x),(-gk));
        nd_x = nd_x + dk;
        k = k + 1;
    end
    fprintf("迭代次数：%d\n",k);
    fprintf("极小点：%s\n",mat2str(double(nd_x)));
    fprintf("极小值：%f\n",fx(nd_x));
end

%% 
function fx = fun(x,n)
fx = (x(1)-5)^2;
for i = 2:n
    fx = fx + (sum(x(1:i))-1)^2;
end
end

function g = gradient1(x,n)
g = zeros(n,1);
g(1) = 2 * (x(1) - 5);
for i = 2:n
    g(1) = g(1) + 2*sum(x(1:i)) - 2;
end
for k = 2:n % 每一个xi
    for i = k:n
        g(k) = g(k) + 2*sum(x(1:i)) - 2;
    end
end
end

function H = hessian1(x,n)
H = zeros(n, n);
for i = 1:n
    for j = 1:n
        H(i, j) = 2*(n-max(i,j)+1);
    end
end
end