%% 修正牛顿法
clear;
close all;
format long;

eps = 1e-6;

for n = [10,100,1000,2000,5000,10000]

    fx = @(x) fun(x,n);
    gx = @(x) gradient1(x,n);
    hx = @(x) hessian1(x,n);
    x0 = 0.01*ones(n,1);

    k = 0;
    delta = 0.5;
    tao = 0.0;
    sigma = 0.4;
    nd_x = x0;

    while 1
        gk = gx(nd_x);
        miuk = norm(gk)^(tao+1);
        if(norm(gk)<=eps)
            break;
        end
        Gk = hx(nd_x);
        dk = -(Gk+miuk*eye(n))\gk;

        m = 0;
        while 1
            if(fx(nd_x+delta^m*dk) <= fx(nd_x)+sigma*delta^m*gk'*dk)
                mk = m;
                break;
            end
            m = m + 1;
        end
        nd_x = nd_x + delta^mk*dk;
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