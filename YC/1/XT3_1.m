clear;
format long;
x = sym('x',[1,2]);
fx = @(x) 3*x(1).^2+2*x(2).^2-4*x(1)-6*x(2);
FX = fx(x);
gfx = gradient(FX); % 梯度
Gfx = hessian(FX); % Hessian矩阵
gfx = matlabFunction(gfx);

% 初始化参数
x0 = [0,3]';
k = 0;
eps = 1e-6;

%% 最速下降法
disp('最速下降法');
ds_x = x0;
while(1)
    gk = gfx(ds_x(1),ds_x(2));
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
fprintf("最小值：%f\n\n", fx(ds_x));

%% 牛顿法
disp('牛顿法');
k = 0;
nd_x = x0;

while(1)
    gk = gfx(nd_x(1),nd_x(2));
    if(norm(gk) <= eps)
        break;
    end
    dk = subs(Gfx,x,num2cell(nd_x'))\(-gk);
    nd_x = nd_x + dk;
    k = k + 1;
end
fprintf("迭代次数：%d\n",k);
fprintf("最小点：%s\n", mat2str(double(nd_x)));
fprintf("最小值：%f\n", fx(nd_x));