clear;
n = 12;
x = sym('x',[1,n]);
fx = @(x) fun(x,n);
gx = @(x) gradient1(x,n);
hx = @(x) hessian1(x,n);

H = hx(ones(1,12))
Hx = matlabFunction(hessian(fx(x)));
Hx(1,1,1,1,1,1,1,1,1,1,1,1)

% gx([1,2,3,4,5,6,7,8,9,10,11,12])'
% Gx = matlabFunction(gradient(fx(x)));
% Gx(1,2,3,4,5,6,7,8,9,10,11,12)

%%
fff = (x(1)-3)^2+(x(1)-3-2*(sum(x(1:2)))^2)^2 +(x(1)-3-2*(sum(x(1:3)))^2)^2;
ffff = diff(fff,x(2));
diff(ffff,x(2))

%%
function fx = fun(x,n)
fx = (x(1)-3)^2;
for i = 2:n
    fx = fx + (x(1) - 3 - 2*sum(x(1:i))^2)^2;
end
end

%%
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

%%
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
                H(i,j) = H(i,j) + 48*sum(x(1:m))^2 - 8*x(1) + 24;
            end
        end
    end
end
end

%%
function He = hessian2(x,n)
for j2=1:n
    f=0;
    for i2=1:n
        if j2==1&&i2==1
            for j3=2:n
                f1=-8*(n-1)*x(1)+26+24*(n-2);%n大于等于2
                f2=16*(sum(x(1:j3))^2)+2*(4*(sum(x(1:j3)))-1)^2;
                if j3==2
                    f=f+f1+f2;
                else
                    f=f+f2;
                end
            end
        else
            if j2==i2
                for j3=n:-1:j2
                    f1=-8*(n-j2+1)*x(1)+24*(n-j2+1);%n大于等于2
                    f2=16*(sum(x(1:j3))^2)+2*(4*(sum(x(1:j3))))^2;
                    if j3==n
                        f=f+f1+f2;
                    else
                        f=f+f2;
                    end
                end
            else
                if j2==1&&i2~=1
                    for j3=n:-1:i2
                        f1=-8*(n-i2+1)*x(1)+24*(n-i2+1);%n大于等于2
                        f2=16*(sum(x(1:j3))^2)+2*(4*sum(x(1:j3))*(4*(sum(x(1:j3)))-1));
                        if j3==n
                            f=f+f1+f2;
                        else
                            f=f+f2;
                        end
                    end
                else
                    if j2~=1&&i2~=1&&i2>j2
                        for j3=n:-1:i2
                            f1=-8*(n-i2+1)*x(1)+24*(n-i2+1);%n大于等于2
                            f2=16*(sum(x(1:j3))^2)+2*(4*sum(x(1:j3)))^2;
                            if j3==n
                                f=f+f1+f2;
                            else
                                f=f+f2;
                            end
                        end
                    end
                end
            end
        end
        He(j2,i2)=f;
        f=0;
    end
end
He=triu(He,1)'+He;
end