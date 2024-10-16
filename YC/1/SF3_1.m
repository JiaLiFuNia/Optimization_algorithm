% 最速下降法

clear;

fx = @(x) 100*(x(1).^2-x(2)).^2+(x(1)-1).^2;
syms x y;
FX = fx([x,y]);
gfx = gradient(FX); % 求梯度
gfx = matlabFunction(gfx);

x = [2,1]';
k = 0;
eps = 1e-6;

while(1)
    gk = gfx(x(1),x(2));
    if(norm(gk) <= eps)
        break;
    end
    dk = -gk;

    beta = 0.5;
    sigma = 0.3;
    m = 0;
    while(m < 20)
        if(fx(x+beta^m*dk) <= fx(x)+sigma*beta^m*gk'*dk)
            mk = m;
            break;
        end
        m = m + 1;
    end
    x = x + beta^mk*dk;
    k = k + 1;
end
k
x
fx(x)