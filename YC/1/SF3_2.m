% 牛顿法

clear;
format long;
x = sym('x',[1,2]);
fx = @(x) 3*x(1).^2+2*x(2).^2-4*x(1)-6*x(2);
FX = fx(x);
gfx = gradient(FX); % 梯度
Gfx = hessian(FX); % Hessian矩阵
gfx = matlabFunction(gfx);

eps = 1e-6;

x0 = [-1,-1]';
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