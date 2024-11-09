%% 牛顿-最速混合法
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
    ds_x = x0;
    while 1
        gk = gx(ds_x);
        if(norm(gk) <= eps)
            break;
        end
        Gk = hx(ds_x);
        dk = Gk\(-gk);
        if(gk'*dk >=0)
            dk = -gk;
        end

        beta = 0.5;
        sigma = 0.3;
        m = 0;
        while 1
            if(fx(ds_x+beta^m*dk) <= fx(ds_x)+sigma*beta^m*gk'*dk)
                mk = m;
                break;
            end
            m = m + 1;
        end
        ds_x = ds_x + beta^mk*dk;
        k = k + 1;
    end
    fprintf("迭代次数：%d\n",k);
    fprintf("极小点：%s\n",mat2str(double(ds_x)));
    fprintf("极小值：%f\n",fx(ds_x));
    
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