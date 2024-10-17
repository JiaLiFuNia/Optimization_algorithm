clear;
format long;

min_value = []; % 极小值
fun_caculate = 0; % 函数值计算次数
gran_caculate = 0; % 梯度计算次数
k = 0; % 迭代次数

ds_res = {}; % 输出结果
nd_res = {};

count = 0;

c = 10;
for n = [10,100,1000,2000,5000,10000]
    count = count + 1;
    x0 = 0.1*ones(n,1); % [0,1]随机数
    fx = @(x) fun(x,c,n);
    gfx = @(x) gradient(x,c,n); % 梯度
    Gfx = @(x) hessian(x,c,n); % Hessian矩阵

    % 初始化参数
    k = 0;
    eps = 1e-4;

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
            m = m + 1;
        end
        ds_x = ds_x + beta^mk*dk;
        fun_caculate = fun_caculate + 2;
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
    min_value = [];
    disp('牛顿法');
    k = 0;
    nd_x = x0;
    tic;
    while(1)
        gk = gfx(nd_x);
        if(norm(gk) <= eps || isnan(nd_x(1)))
            break;
        end
        dk = lsqminnorm(Gfx(nd_x),(-gk));
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
figure(1)
set(gcf,'Position',[100,200,800,450]);

subplot(1,2,1);
for i = 1:size(ds_res,1)
    plot(log10(ds_res{i,4}),'LineWidth',1.5);
    hold on;
end
title('函数值随迭代次数的变化(最速下降法)')
xlabel('迭代次数');
ylabel('$log_{10}(f(x_{min}))$','Interpreter','latex');
% [10,100,1000,2000,5000,10000]
h = legend({'$n=10$','$n=100$','$n=1000$','$n=2000$','$n=5000$','$n=10000$'});
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
h = legend({'$n=10$','$n=100$','$n=1000$','$n=2000$','$n=5000$','$n=10000$'});
set(h,'Interpreter','latex');
set(gca,'Position',[0.56,0.1,0.42,0.85]);
hold off;

% figure(2)
figure(2)
plot(log10([ds_res{:,3}]),'LineWidth',1.5);
hold on;
plot(log10([nd_res{:,3}]),'LineWidth',1.5);
xticks(1:1:4) % x的范围
xticklabels({'10','100','1000','2000','5000','10000'}) % 修改x的标签
title('运行时间随n的变化')
legend({'最速下降法','牛顿法'});
xlabel('$\bf{n}$','Interpreter','latex');
ylabel('$\bf{log_{10}(t)}$','Interpreter','latex');

%% fx
function fx = fun(x,c,n)
fx = 0;
for i = 1:n
    fx = fx + exp(x(i))-x(i);
end
end

%% 梯度
function grad = gradient(x,c,n)
grad = zeros(n,1);
for i = 1:n
    grad(i) = exp(x(i)) - 1;
end
end

%% Hessian
function H = hessian(x,c,n)
H = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i == j
            H(i,j) = exp(x(i));
        else
            H(i,j) = 0;
        end
    end
end
end