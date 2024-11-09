%% FR共轭梯度法  
clear;
close all;
format long;

eps = 1e-6;
delta = 0.5;
sigma = 0.4;

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
fx = (x(1)-3)^2;
for i = 2:n
    fx = fx + (x(1) - 3 - 2*sum(x(1:i))^2)^2;
end
end

function g = gradient1(x,n)
g = zeros(n,1);
g(1) = 2 * (x(1) - 3);
for i = 2:n
    temp = 2 * sum(x(1:i))^2 - x(1) + 3;
    g(1) = g(1) + 2 * temp * (4*sum(x(1:i))-1);
end
for k = 2:n % 每一个xi
    for i = k:n
        temp = 2 * sum(x(1:i))^2 - x(1) + 3;
        g(k) = g(k) + 2 * temp * 4 * sum(x(1:i));
    end
end
end

function H = hessian1(x,n)
H = zeros(n, n);
H(1,1) = 2;
for i = 2:n
    H(1,1) = H(1,1) + 16*sum(x(1:i))^2 + 2*(4*sum(x(1:i))-1)^2 - 8*x(1) + 24;
end
for i = 1:n
    for j = 1:n
        if i == 1 && j == 1
            continue;
        end
        for m = max(i,j):n
            if i == 1 || j == 1 % 第1行和第1列并排除H(1,1)
                H(i,j) = H(i,j) + 16*sum(x(1:m))^2 + 8*sum(x(1:m))*(4*sum(x(1:m))-1) - 8*x(1) + 24;
            else
                H(i,j) = H(i,j) + 48*sum(x(1:m))^2 - 8*x(1) + 24;
            end
        end
    end
end
end