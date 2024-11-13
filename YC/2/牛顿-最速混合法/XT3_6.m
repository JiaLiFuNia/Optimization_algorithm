%% 牛顿-最速混合法
clear;
close all;
format long;

load('k3_6.mat');
load('time3_6.mat')

eps = 1e-6;

i = 0;
for n = [10,100,1000,2000,5000,10000]
    tic;
    x0 = 0.9*ones(n,1);
    for c = 10
        [n,c]
        fx = @(x) fun(x,c,n);
        gx = @(x) gradient(x,c,n); % 梯度
        hx = @(x) hessian(x,c,n); % Hessian矩阵

        k = 0;
        ds_x = x0;
        while k <= 20000 % 避免某些情况下迭代时间太长
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
    i = i + 1;
    K(2,i) = k;
    TIME(2,i) = toc;
end
save('k3_6.mat',"K");
save('time3_6.mat',"TIME");

%% fx
function fx = fun(x,c,n)
fx = 0;
for i = 1:n-1
    fx = fx + c*(x(i)^2-x(i+1)).^2 + (x(i)-1).^2;
end
end

%% 梯度
function grad = gradient(x,c,n)
grad = zeros(n,1);
for i = 1:n
    if i < n
        grad(i) = 4*c*x(i)*(x(i)^2 - x(i+1)) + 2*(x(i) - 1);
    end
    if i > 1
        grad(i) = grad(i) - 2*c*(x(i-1)^2 - x(i));
    end
end
end

%% Hessian
function H = hessian(x,c,n)
H = zeros(n,n);
for i = 1:n
    if i < n
        H(i,i) = 12*c*x(i)^2 - 4*c*x(i+1) + 2; % 二阶导
        H(i,i+1) = -4*c*x(i);                   % 对 x(i) 和 x(i+1) 的混合导数
    end
    if i > 1
        H(i,i) = H(i,i) + 2*c;
        H(i,i-1) = H(i-1,i);  % 对 x(i) 和 x(i-1) 的混合导数
    end
end
end