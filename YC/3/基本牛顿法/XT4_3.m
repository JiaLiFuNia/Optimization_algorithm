%% 基本牛顿法
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

    i = i + 1;
    K(1,i) = k;
    TIME(1,i) = toc;
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