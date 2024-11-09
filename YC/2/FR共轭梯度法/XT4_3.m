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
fx = 0;
for i = 1:n-1
    fx = fx+cos(-0.5*x(i+1)+x(i)^2);
end
end

function grad = gradient1(x,n)
grad = zeros(n,1);
for i = 1:n
    if i < n
        grad(i) =(-2)*x(i)*sin(x(i)^2-0.5*x(i+1));
    end
    if i > 1
        grad(i) = grad(i)+(0.5)*sin(x(i-1)^2-0.5*x(i));
    end
end
end

function H = hessian1(x,n)
H = zeros(n,n);
for i = 1:n
    if i < n
        H(i,i) = (-2)*sin(x(i)^2-0.5*x(i+1))-4*x(i)^2*cos(x(i)^2-0.5*x(i+1));
        H(i,i+1) = x(i)*cos(x(i)^2-0.5*x(i+1));                 
    end
    if i > 1
        H(i,i) = H(i,i)-0.25*cos(x(i-1)^2-0.5*x(i));
        H(i,i-1) = H(i-1,i);  
    end
end
end