%% 牛顿-最速混合法
clear;
close all;
format long;

eps = 1e-6;

for n = [10,100,1000,2000,5000,10000]

    fx = @(x) fun(x,n);
    gx = @(x) gradient1(x,n);
    hx = @(x) hessian1(x,n);
    x0 = 2*ones(n,1);

    k = 0;
    ds_x = x0;
    while 1
        gk = gx(ds_x);
        if(norm(gk) <= eps)
            break;
        end
        Gk = hx(ds_x);
        dk = Gk\(-gk);
        if(gk'*dk >=0)
            dk = -gk;
        end

        beta = 0.5;
        sigma = 0.3;
        m = 0;
        while 1
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
    fprintf("极小点：%s\n",mat2str(double(ds_x)));
    fprintf("极小值：%f\n",fx(ds_x));
    
end

%% 
function fx = fun(x,n)
fx = 0;
for i = 1:n-1
    fx = fx+sin(-0.5*x(i+1)+x(i)^2);
end
end    

function grad = gradient1(x,n)
grad = zeros(n,1);
for i = 1:n
    if i < n
        grad(i) =(2)*x(i)*cos(x(i)^2-0.5*x(i+1));
    end
    if i > 1
        grad(i) = grad(i)+(-0.5)*cos(x(i-1)^2-0.5*x(i));
    end
end
end

function H = hessian1(x,n)
H = zeros(n,n);
for i = 1:n
    if i < n
        H(i,i) = (2)*cos(x(i)^2-0.5*x(i+1))-4*x(i)^2*sin(x(i)^2-0.5*x(i+1));
        H(i,i+1) = x(i)*sin (x(i)^2-0.5*x(i+1));                 
    end
    if i > 1
        H(i,i) = H(i,i)-0.25*sin(x(i-1)^2-0.5*x(i));
        H(i,i-1) = H(i-1,i);  
    end
end
end