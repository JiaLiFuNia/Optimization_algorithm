clear;
n = 12;
x = sym('x',[1,n]);
fx = @(x) fun(x,n);
gx = @(x) gradient1(x,n);
hx = @(x) hessian1(x,n);

H = hx([1,2,3,4,5,6,7,8,9,10,11,12])
Hx = matlabFunction(hessian(fx(x)));
Hx(1,2,3,4,5,6,7,8,9,10,11,12)

% G = gradient(fx(x))
%
% gx([1,2,3,4,5,6,7,8,9,10,11,12])
% Gx = matlabFunction(G);
% Gx(1,2,3,4,5,6,7,8,9,10,11,12)
%
%%
fff = (x(1)-3-2*(sum(x(1:2)))^2)^2 + (x(1)-3-2*(sum(x(1:3)))^2)^2 + (x(1)-3-2*(sum(x(1:4)))^2)^2 + (x(1)-3-2*(sum(x(1:5)))^2)^2;
ffff = diff(fff,x(3));
fffff = diff(ffff,x(2))

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
    temp = x(1) - 3 - 2 * sum(x(1:i))^2;
    g(1) = g(1) + (-2) * temp * (4*sum(x(1:i))-1);
end
for k = 2:n
    for i = k:n
        temp = x(1) - 3 - 2 * sum(x(1:i))^2;
        g(k) = g(k) + (-2) * temp * 4 * sum(x(1:i));
    end
end
end

function H = hessian1(x,n)
H = zeros(n, n);
H(1,1) = 2;
for i = 2:n
    H(1,1) = H(1,1) + 16*sum(1:i)^2 + 2*(4*sum(x(1:i))-1)^2 - 8*x(1) + 24;
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
                H(i,j) = H(i,j) + 32*sum(x(1:m))^2 - 8*x(1) + 24;
            end
        end
    end
end
end