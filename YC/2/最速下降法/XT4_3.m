%% 最速下降法
clear;
close all;
format long;

load('k4_3.mat');
load('time4_3.mat')

eps = 1e-4;

i = 0;
for n = [10,100,1000,2000]

    tic;
    x0 = 2*ones(n,1);
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
   
    i = i + 1;
    K(3,i) = k;
    TIME(3,i) = toc;
    fprintf("迭代时间：%f\n\n", TIME(2,i));
end
save('k4_3.mat',"K");
save('time4_3.mat',"TIME");

%% 
function fx = fun(x,n)
fx = 0;
for i = 1:n-1
    fx = fx+cos(-0.5*x(i+1)+x(i)^2);
end
end

function grad = gradient1(x,n)
grad = zeros(n,1);
for i = 1:n
    if i < n
        grad(i) =(-2)*x(i)*sin(x(i)^2-0.5*x(i+1));
    end
    if i > 1
        grad(i) = grad(i)+(0.5)*sin(x(i-1)^2-0.5*x(i));
    end
end
end

function H = hessian1(x,n)
H = zeros(n,n);
for i = 1:n
    if i < n
        H(i,i) = (-2)*sin(x(i)^2-0.5*x(i+1))-4*x(i)^2*cos(x(i)^2-0.5*x(i+1));
        H(i,i+1) = x(i)*cos(x(i)^2-0.5*x(i+1));                 
    end
    if i > 1
        H(i,i) = H(i,i)-0.25*cos(x(i-1)^2-0.5*x(i));
        H(i,i-1) = H(i-1,i);  
    end
end
end