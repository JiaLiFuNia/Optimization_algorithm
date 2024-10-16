% 牛顿法

clear;
fx = @(x) 3*x(1).^2+2*x(2).^2-4*x(1)-6*x(2);
syms x y;
FX = fx([x,y]);
gfx = gradient(FX); % 求梯度
Gfx = hessian(FX); % Hessian矩阵

if isempty(symvar(Gfx))
    flag = 1;
else
    flag = 0;
end

gfx = matlabFunction(gfx);
Gfx = matlabFunction(Gfx);

x = [-1,-1]';
k = 0;
eps = 1e-6;

while(1)
    gk = gfx(x(1),x(2));
    if(norm(gk) <= eps)
        break;
    end
    if flag == 1
        dk = Gfx()\(-gk);
    else
        dk = Gfx(x(1),x(2))\(-gk);
    end
    x = x + dk;
    k = k + 1;
end
k
x
fx(x)