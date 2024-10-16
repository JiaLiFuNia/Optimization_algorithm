clear;
close all;

format long;
% 定义函数以及初始点
phi = @(x) x.^2-x-1;
fplot(phi);
a0 = -1;b0 = 1;
xlim([a0,b0]);
t = (sqrt(5)-1)/2;
eps = 1e-5;
delta = 0.05;

a = [];
b = [];
p = [];
q = [];

% 初始化向量
a(1) = a0;
b(1) = b0;
p(1) = a(1) + (1-t)*(b(1)-a(1));
q(1) = a(1) + t*(b(1)-a(1));
h = abs(a(1) - b(1));

i = 0;
phip = phi(p(1));
phiq = phi(q(1));

while(abs(q(i+1)-a(i+1)) > eps || abs(b(i+1)-p(i+1)) > eps || h > delta)
    if(phip <= phiq)
        a(i+2) = a(i+1);
        b(i+2) = q(i+1);
        phiq = phip;
        q(i+2) = p(i+1);
        p(i+2) = a(i+2) + (1-t)*(b(i+2)-a(i+2));
        phip = phi(p(i+2));
    else
        a(i+2) = p(i+1);
        b(i+2) = b(i+1);
        phip = phiq;
        p(i+2) = q(i+1);
        q(i+2) = a(i+2) + t*(b(i+2)-a(i+2));
        phiq = phi(q(i+2));
    end
    h = abs(a(i+2)-b(i+2));
    i = i + 1;
end
if(phip <= phiq)
    minNode = p(i+1);
    min = phip;
else
    minNode = q(i+1);
    min = phiq;
end
disp([minNode,min,a(i),b(i)]);