%% 修正牛顿法
clear;
close all;
format long;

% x = sym('x',[1,n]);

% H = hx(ones(1,n))
% Hx = matlabFunction(hessian(fx(x)));
% Hx(1,1,1)

% gx([1,2,3,4,5,6,7,8,9,10,11,12])'
% Gx = matlabFunction(gradient(fx(x)));
% Gx(1,2,3,4,5,6,7,8,9,10,11,12)

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
% fff = (x(1)-3)^2+(x(1)-3-2*(sum(x(1:2)))^2)^2;
% ffff = diff(fff,x(1));
% fffff = diff(ffff,x(1))
% fc = matlabFunction(fffff)

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