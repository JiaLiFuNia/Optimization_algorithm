% 对称秩1算法
clear;
format long;

eps = 1e-6;
delta = 0.55;
sigma = 0.4;
n = 2;

H = eye(n);
gx = @(x) gradient(x);
fx = @(x) fun(x);

x0 = [0.5;0.5];
k = 0;
gk = [];
dk = [];

dcz_x = [];
dcz_x(:,1) = x0;
while 1
    gk(:,k+1) = gx(dcz_x(:,k+1));
    if(norm(gk(:,k+1)) <= eps)
        break;
    end
    H;
    dk(:,k+1) = -H*gk(:,k+1);

    m = 0;
    while m <= 20
        if(fx(dcz_x(:,k+1)+delta^m*dk(:,k+1)) <= fx(dcz_x(:,k+1))+sigma*delta^m*gk(:,k+1)'*dk(:,k+1))
            mk = m;
            break;
        else
            m = m + 1;
        end
    end
    dcz_x(:,k+2) = dcz_x(:,k+1) + delta^mk*dk(:,k+1);
    sk = delta^mk*dk(:,k+1);
    yk = gx(dcz_x(:,k+2)) - gk(:,k+1);
    H = H + (sk-H*yk)*(sk-H*yk)'/((sk-H*yk)'*yk);
    k = k + 1;
end
fprintf("迭代次数：%d\n",k);
fprintf("极小点：%s\n",mat2str(dcz_x(:,end)));
fprintf("极小值：%f\n",fx(dcz_x(:,end)));


%%
% syms x [1 2]
% f = 100*(x(1)^2-x(2))^2+(x(1)-1)^2;
% ff = diff(f,x(2))


%%
function f = fun(x)
f = 100*(x(1)^2-x(2))^2+(x(1)-1)^2;
end

function g = gradient(x)
g = [2*x(1) - 400*x(1)*(- x(1)^2 + x(2)) - 2;- 200*x(1)^2 + 200*x(2)];
end
