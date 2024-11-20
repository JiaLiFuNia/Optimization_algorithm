%% BFGS
clear;
format long;

eps = 1e-6;
delta = 0.55;
sigma = 0.4;


for n = [10,100,1000,2000,5000,10000,100000]
    k = 0;

    gx = @(x) gradient1(x,n);
    fx = @(x) fun(x,n);

    B = eye(n);
    x0 = zeros(n,1);
    x0(1:2) = [2.9,-2.9];

    gk = [];
    dk = [];
    dcz_x = [];
    dcz_x(:,1) = x0;

    while 1
        gk(:,k+1) = gx(dcz_x(:,k+1));
        if(norm(gk(:,k+1)) <= eps)
            break;
        end
        dk(:,k+1) = -B\gk(:,k+1);

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
        if yk'*sk > 0
            B = B - B*(sk*sk')*B/(sk'*B*sk) + (yk*yk')/(yk'*sk);
        end
        k = k + 1;
    end
    fprintf("n:%d\n",n);
    fprintf("迭代次数：%d\n",k);
    fprintf("极小点：%s\n",mat2str(dcz_x(:,end)));
    fprintf("极小值：%f\n\n",fx(dcz_x(:,end)));
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
