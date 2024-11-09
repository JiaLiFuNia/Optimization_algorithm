%% FR共轭梯度法
clear;
close all;
format long;

eps = 1e-6;

for n = [10,100,1000,2000,5000,10000]
    fx = @(x) fun(x,n);
    gx = @(x) gradient1(x,n);
    hx = @(x) hessian1(x,n);
    x0 = 0.01*ones(n,1);
  
    k = 0;
    ge_x = x0;
    gk = [];
    dk = [];

    while 1
        gk(:,k+1) = gx(ge_x);
        if(norm(gk(:,k+1)) <= eps)
            break;
        end
        if k == 0
            dk(:,k+1) = -gk(:,k+1);
        else
            beta_k = gk(:,k+1)'*gk(:,k+1)/(gk(:,k)'*gk(:,k));
            dk(:,k+1) = -gk(:,k+1) + beta_k*dk(:,k);
        end

        m = 0;
        while 1
            if(fx(ge_x+delta^m*dk(:,k+1)) <= fx(ge_x)+sigma*delta^m*gk(:,k+1)'*dk(:,k+1))
                mk = m;
                break;
            end
            m = m + 1;
        end

        ge_x = ge_x + delta^mk*dk(:,k+1);
        k = k + 1;
    end
    fprintf("迭代次数：%d\n",k);
    fprintf("极小点：%s\n",mat2str(double(ge_x)));
    fprintf("极小值：%f\n",fx(ge_x));

end

%%
function fx = fun(x,n)
fx = (x(1)-5)^2;
for i = 2:n
    fx = fx + (sum(x(1:i))-1)^2;
end
end

function g = gradient1(x,n)
g = zeros(n,1);
g(1) = 2 * (x(1) - 5);
for i = 2:n
    g(1) = g(1) + 2*sum(x(1:i)) - 2;
end
for k = 2:n % 每一个xi
    for i = k:n
        g(k) = g(k) + 2*sum(x(1:i)) - 2;
    end
end
end

function H = hessian1(x,n)
H = zeros(n, n);
for i = 1:n
    for j = 1:n
        H(i, j) = 2*(n-max(i,j)+1);
    end
end
end