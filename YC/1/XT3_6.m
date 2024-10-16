clear;
format long;

for n = [10,100,1000,2000,5000,10000]
    x0 = 0.9*ones(n,1); % [0,1]随机数
    for c = [1,10,100,1000]
        n,c
        fx = @(x) fun(x,c,n);
        gfx = @(x) gradient(x,c,n); % 梯度
        Gfx = @(x) hessian(x,c,n); % Hessian矩阵

        % 初始化参数
        k = 0;
        eps = 1e-4;

        %% 最速下降法
        disp('最速下降法');
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
            while(m < 20)
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
        fx(ds_x)
        fprintf("最小值：%f\n\n", fx(ds_x));

        %% 牛顿法
        disp('牛顿法');
        k = 0;
        nd_x = x0;

        while(1)
            gk = gfx(nd_x);
            if(norm(gk) <= eps || isnan(nd_x(1)))
                break;
            end
            dk = lsqminnorm(Gfx(nd_x),(-gk));
            nd_x = nd_x + dk;
            k = k + 1;
        end
        fprintf("迭代次数：%d\n",k);
        fprintf("最小点：%s\n", mat2str(double(nd_x)));
        fx(ds_x)
        fprintf("最小值：%f\n", fx(nd_x));
    end
end

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