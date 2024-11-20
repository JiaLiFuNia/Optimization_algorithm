%% 修正牛顿法
clear;
close all;
format long;

load('k4_4.mat');
load('time4_4.mat')

eps = 1e-4;

i = 0;
for n = [10,100,1000,2000]
    for x = 2.40:0.01:2.49

    tic;

    fx = @(x) fun(x,n);
    gx = @(x) gradient1(x,n);
    hx = @(x) hessian1(x,n);
    x0 = x*ones(n,1);

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

    i = i + 1;
    K(5,i) = k;
    TIME(5,i) = toc;
    fprintf("迭代时间：%f\n\n", TIME(2,i));
    end
    i = i + 1;
end
save('k4_4.mat',"K");
save('time4_4.mat',"TIME");

%% 
function fx = fun(x,n)
fx = 0;
for i = 1:n-1
    fx = fx+sin(-0.5*x(i+1)+x(i)^2);
end
end    

function grad = gradient1(x,n)
grad = zeros(n,1);
for i = 1:n
    if i < n
        grad(i) =(2)*x(i)*cos(x(i)^2-0.5*x(i+1));
    end
    if i > 1
        grad(i) = grad(i)+(-0.5)*cos(x(i-1)^2-0.5*x(i));
    end
end
end

function H = hessian1(x,n)
H = zeros(n,n);
for i = 1:n
    if i < n
        H(i,i) = (2)*cos(x(i)^2-0.5*x(i+1))-4*x(i)^2*sin(x(i)^2-0.5*x(i+1));
        H(i,i+1) = x(i)*sin (x(i)^2-0.5*x(i+1));                 
    end
    if i > 1
        H(i,i) = H(i,i)-0.25*sin(x(i-1)^2-0.5*x(i));
        H(i,i-1) = H(i-1,i);  
    end
end
end