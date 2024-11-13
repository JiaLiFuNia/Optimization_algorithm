clear;
close all;
n = 10;
x = sym('x',[1,n]);
fx = @(x) fun(x,n);
gx = @(x) gradient1(x,n);
hx = @(x) hessian1(x,n);

H = hx([1,2,3,4,5,6,7,8,9,10])
Hx = matlabFunction(hessian(fx(x)));
h = Hx(1,2,3,4,5,6,7,8,9,10)

% gx([1,2,3,4,5,6,7,8,9,10,11,12])
% Gx = matlabFunction(gradient(fx(x)));
% Gx(1,2,3,4,5,6,7,8,9,10,11,12)

%%
x = sym('x',[1,n]);
fff = cos(-0.5*x(2)+x(1)^2);
ffff = diff(fff,x(1))
fffff = diff(ffff,x(1))
% fc = matlabFunction(fffff)

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
        grad(i) = (-2)*x(i)*sin(x(i)^2-0.5*x(i+1));
    end
    if i > 1
        grad(i) = grad(i) + (0.5)*sin(x(i-1)^2-0.5*x(i));
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
        H(i,i) = H(i,i) - 0.25*cos(x(i-1)^2-0.5*x(i));
        H(i,i-1) = H(i-1,i);  
    end
end
end