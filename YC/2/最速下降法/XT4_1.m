%% 最速下降法
clear;
close all;
format long;

eps = 1e-4;

for n = [10,100,1000,2000,5000,10000]

    x0 = 0.1*ones(n,1);
    fx = @(x) fun(x,n);
    gfx = @(x) gradient1(x,n); % 梯度
    Gfx = @(x) hessian1(x,n); % Hessian矩阵

    k = 0;
    ds_x = x0;
    while(1)
        gk = gfx(ds_x);
        if(norm(gk) <= eps)
            break;
        end
        dk = -gk;

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
    fprintf("最小点：%s\n", mat2str(ds_x));
    fprintf("最小值：%f\n\n", fx(ds_x));
   
end

%% 
function fx = fun(x,n)
fx = (x(1)-3)^2;
for i = 2:n
    fx = fx + (x(1) - 3 - 2*sum(x(1:i))^2)^2;
end
end

function g = gradient1(x,n)
g = zeros(n,1);
g(1) = 2 * (x(1) - 3);
for i = 2:n
    temp = 2 * sum(x(1:i))^2 - x(1) + 3;
    g(1) = g(1) + 2 * temp * (4*sum(x(1:i))-1);
end
for k = 2:n % 每一个xi
    for i = k:n
        temp = 2 * sum(x(1:i))^2 - x(1) + 3;
        g(k) = g(k) + 2 * temp * 4 * sum(x(1:i));
    end
end
end

function H = hessian1(x,n)
H = zeros(n, n);
H(1,1) = 2;
for i = 2:n
    H(1,1) = H(1,1) + 16*sum(x(1:i))^2 + 2*(4*sum(x(1:i))-1)^2 - 8*x(1) + 24;
end
for i = 1:n
    for j = 1:n
        if i == 1 && j == 1
            continue;
        end
        for m = max(i,j):n
            if i == 1 || j == 1 % 第1行和第1列并排除H(1,1)
                H(i,j) = H(i,j) + 16*sum(x(1:m))^2 + 8*sum(x(1:m))*(4*sum(x(1:m))-1) - 8*x(1) + 24;
            else
                H(i,j) = H(i,j) + 48*sum(x(1:m))^2 - 8*x(1) + 24;
            end
        end
    end
end
end