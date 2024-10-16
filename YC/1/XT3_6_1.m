clear;
close all;
format long;

min_value = []; % 极小值
fun_caculate = 0; % 函数值计算次数
gran_caculate = 0; % 梯度计算次数
k = 0; % 迭代次数

res = {}; % 输出结果

%% 初始化参数
eps = 1e-4; % 精度
n = 10000;
c = 10;
x0 = 0.9*ones(n,1); % 初始值

fx = @(x) fun(x,c,n);
gfx = @(x) gradient(x,c,n); % 梯度
Gfx = @(x) hessian(x,c,n); % Hessian矩阵

%% 最速下降法
disp('最速下降法');
ds_x = x0;
tic;
while(1)
    gk = gfx(ds_x);
    gran_caculate = gran_caculate + 1;
    if(norm(gk) <= eps)
        break;
    end
    dk = -gk;

    beta = 0.5;
    sigma = 0.3;
    m = 0;
    while(m < 20)
        if(fx(ds_x+beta^m*dk) <= fx(ds_x)+sigma*beta^m*gk'*dk)
            mk = m;
            break;
        end
        fun_caculate = fun_caculate + 2;
        m = m + 1;
    end
    ds_x = ds_x + beta^mk*dk;
    k = k + 1;
    min_value(k) = fx(ds_x);
end
running_time = toc;
fprintf("迭代次数：%d\n",k);
fprintf("最小点：%s\n", mat2str(ds_x));
fprintf("最小值：%f\n\n", fx(ds_x));

res{1,1} = k;
res{1,2} = fx(ds_x);
res{1,3} = running_time;
res{1,4} = min_value;
subplot(1,2,1);
plot(log10(min_value),'LineWidth',1.5);
title('函数值随迭代次数的变化(最速下降法)')
xlabel('迭代次数');
ylabel('$log_{10}(f(x_{min}))$','Interpreter','latex');
set(gca,'Position',[0.07,0.1,0.42,0.85])

%% 牛顿法
disp('牛顿法');
k = 0;
nd_x = x0;
min_value = [];
tic;
while(1)
    gk = gfx(nd_x);
    if(norm(gk) <= eps)
        break;
    end
    dk = lsqminnorm(Gfx(nd_x),(-gk));
    nd_x = nd_x + dk;
    k = k + 1;
    min_value(k) = fx(nd_x);
end
running_time = toc;
fprintf("迭代次数：%d\n",k);
fprintf("最小点：%s\n", mat2str(double(nd_x)));
fprintf("最小值：%f\n", fx(nd_x));

res{2,1} = k;
res{2,2} = fx(nd_x);
res{2,3} = running_time;
res{2,4} = min_value;
subplot(1,2,2);
plot(log10(min_value),'LineWidth',1.5);
title('函数值随迭代次数的变化(牛顿法)')
xlabel('迭代次数');
ylabel('$log_{10}(f(x_{min}))$','Interpreter','latex');
set(gca,'Position',[0.56,0.1,0.42,0.85])
set(gcf,'Position',[100,200,800,450]);

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