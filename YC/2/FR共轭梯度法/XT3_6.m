%% FR共轭梯度法
clear;
close all;
format long;

load('k3_6.mat');
load('time3_6.mat')

eps = 1e-6;
delta = 0.5;
sigma = 0.4;

i = 0;
for n = [10,100,1000,2000,5000,10000]
    tic;
    x0 = 0.9*ones(n,1);
    for c = 10
        [n,c]
        fx = @(x) fun(x,c,n);
        gx = @(x) gradient(x,c,n); % 梯度
        hx = @(x) hessian(x,c,n); % Hessian矩阵

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
    i = i + 1;
    K(4,i) = k;
    TIME(4,i) = toc;
end
save('k3_6.mat',"K");
save('time3_6.mat',"TIME");

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