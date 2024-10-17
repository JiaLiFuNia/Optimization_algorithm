clear;
format long;
x = sym('x',[1,2]);
fx = @(x) 100*(x(2)-x(1).^2).^2+(x(1)-1)^2;
FX = fx(x);
gfx = gradient(FX); % 梯度
Gfx = hessian(FX); % Hessian矩阵
gfx = matlabFunction(gfx);

min_value = []; % 极小值
fun_caculate = 0; % 函数值计算次数
gran_caculate = 0; % 梯度计算次数
ds_res = {}; % 输出结果
nd_res = {};
count = 0;

X = [[-1:0.1:-0.1]',[1:0.1:1.9]'];
for count = 1:size(X,1)
    % 初始化参数
    x0 = X(count,:)';
    k = 0;
    eps = 1e-6;

    %% 最速下降法
    disp('最速下降法');
    ds_x = x0;
    tic;
    while(1)
        gk = gfx(ds_x(1),ds_x(2));
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
            fun_caculate = fun_caculate + 2;
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
    ds_res{count,7} = x0;

    %% 牛顿法
    disp('牛顿法');
    k = 0;
    nd_x = x0;
    tic;
    fun_caculate = 0; % 函数值计算次数
    gran_caculate = 0; % 梯度计算次数
    min_value = [];
    while(1)
        gk = gfx(nd_x(1),nd_x(2));
        if(norm(gk) <= eps)
            break;
        end
        dk = subs(Gfx,x,nd_x')\(-gk);
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
    nd_res{count,7} = x0;

end

%% figure(1)
close all;
set(gcf,'Position',[100,200,800,450]);

subplot(1,2,1);
for i = 1:size(ds_res,1)
    plot(log10(ds_res{i,4}),'LineWidth',1.5);
    hold on;
end
title('不同初值点函数值随迭代次数的变化(最速下降法)')
xlabel('迭代次数');
ylabel('$log_{10}(f(x))$','Interpreter','latex');
legend({'(-1,1)','(-0.9,1.1)','(-0.8,1.2)','(-0.7,1.3)','(-0.6,1.4)','(-0.5,1.5)','(-0.4,1.6)','(-0.3,1.7)','(-0.2,1.8)','(-0.1,1.9)'});
set(gca,'Position',[0.07,0.1,0.42,0.85]);
hold off;

subplot(1,2,2);
for i = 1:size(nd_res,1)
    plot(log10(nd_res{i,4}),'LineWidth',1.5);
    hold on;
end
title('不同初值点函数值随迭代次数的变化(牛顿法)')
xlabel('迭代次数');
ylabel('$log_{10}(f(x))$','Interpreter','latex');
legend({'(-1,1)','(-0.9,1.1)','(-0.8,1.2)','(-0.7,1.3)','(-0.6,1.4)','(-0.5,1.5)','(-0.4,1.6)','(-0.3,1.7)','(-0.2,1.8)','(-0.1,1.9)'});
set(gca,'Position',[0.56,0.1,0.42,0.85]);
hold off;

% figure(2)
figure(2)
N = [10,100,1000,2000];
plot([ds_res{:,3}],'LineWidth',1.5);
hold on;
plot([nd_res{:,3}],'LineWidth',1.5);
xticks(1:1:10) % x的范围
xticklabels({'(-1,1)','(-0.9,1.1)','(-0.8,1.2)','(-0.7,1.3)','(-0.6,1.4)','(-0.5,1.5)','(-0.4,1.6)','(-0.3,1.7)','(-0.2,1.8)','(-0.1,1.9)'}) % 修改x的标签
title('运行时间随初值点的变化')
legend({'最速下降法','牛顿法'});
xlabel('$\bf{n}$','Interpreter','latex');
ylabel('$\bf{t}$','Interpreter','latex');

%% 
for i = 1:size(X,1)
    disp(['(',num2str(X(i,1)),',',num2str(X(i,2)),')']);
end