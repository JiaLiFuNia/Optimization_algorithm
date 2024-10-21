clear;
close all;
format long;

min_value = []; % 极小值
fun_caculate = 0; % 函数值计算次数
gran_caculate = 0; % 梯度计算次数
k = 0; % 迭代次数

ds_res = {}; % 输出结果
nd_res = {};

%% 初始化参数
eps = 1e-4; % 精度
count = 0;
for n = [10,100,1000,2000]
    count = count + 1;
    c = 10;
    x0 = 0.85*ones(n,1); % 初始值

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

    ds_res{count,1} = k;
    ds_res{count,2} = fx(ds_x);
    ds_res{count,3} = running_time;
    ds_res{count,4} = min_value;
    ds_res{count,5} = fun_caculate;
    ds_res{count,6} = gran_caculate;

    %% 牛顿法
    fun_caculate = 0; % 函数值计算次数
    gran_caculate = 0; % 梯度计算次数
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
        tic
        dk = lsqminnorm(Gfx(nd_x),(-gk));
        toc
        % dk = Gfx(nd_x)\(-gk);
        nd_x = nd_x + dk;
        k = k + 1;
        min_value(k) = fx(nd_x);
        gran_caculate = gran_caculate + 1;
    end
    running_time = toc;
    fprintf("迭代次数：%d\n",k);
    fprintf("最小点：%s\n", mat2str(double(nd_x)));
    fprintf("最小值：%f\n", fx(nd_x));

    nd_res{count,1} = k;
    nd_res{count,2} = fx(nd_x);
    nd_res{count,3} = running_time;
    nd_res{count,4} = min_value;
    nd_res{count,5} = fun_caculate;
    nd_res{count,6} = gran_caculate;

end

%% figure(1)
f = figure(1);
set(gcf,'Position',[100,200,800,450]);

subplot(1,2,1);
for i = 1:size(ds_res,1)
    plot(log10(ds_res{i,4}),'LineWidth',1.5);
    hold on;
end
title('函数值随迭代次数的变化(最速下降法)')
xlabel('迭代次数');
ylabel('$log_{10}(f(x_{min}))$','Interpreter','latex');
h = legend({'$n=10$','$n=100$','$n=1000$','$n=2000$'});
set(h,'Interpreter','latex');
set(gca,'Position',[0.07,0.1,0.42,0.85]);
hold off;

subplot(1,2,2);
for i = 1:size(nd_res,1)
    plot(log10(nd_res{i,4}),'LineWidth',1.5);
    hold on;
end
title('函数值随迭代次数的变化(牛顿法)')
xlabel('迭代次数');
ylabel('$log_{10}(f(x))$','Interpreter','latex');
h = legend({'$n=10$','$n=100$','$n=1000$','$n=2000$'});
set(h,'Interpreter','latex');
set(gca,'Position',[0.56,0.1,0.42,0.85]);
hold off;
saveas(f, '1', 'svg');

%% figure(2)
figure(2)
N = [10,100,1000,2000];
plot(log10([ds_res{:,3}]),'LineWidth',1.5);
hold on;
plot(log10([nd_res{:,3}]),'LineWidth',1.5);
xticks(1:1:4) % x的范围
xticklabels({'10','100','1000','2000'}) % 修改x的标签
title('运行时间随n的变化')
legend({'最速下降法','牛顿法'});
xlabel('$\bf{n}$','Interpreter','latex');
ylabel('$\bf{log_{10}(t)}$','Interpreter','latex');

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