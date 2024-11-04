f = @(u) 1/2 * u;
u0 = 1;
h = 0.1;
real_u = @(t) exp(1/2 * t);
% 真实解
rt = linspace(0, 1, 101);
ru = real_u(rt);
% Euler法
t1 = linspace(0, 1, 1/h);
u1 = zeros(size(t1));
u1(1) = u0;
for i = 2:length(t1)
    u1(i) = u1(i-1) + h * f(u1(i-1));
end
% 梯形法（方式一）
t2 = linspace(0, 1, 1/h);
u2 = zeros(size(t2));
u2(1) = u0;
for i = 2:length(t2)
    u2(i) = u2(i-1) + h * f(u2(i-1));
    u2(i) = u2(i-1) + h/2 * (f(u2(i-1)) + f(u2(i)));
end
% 梯形法（方式二）
t3 = linspace(0, 1, 1/h);
u3 = zeros(size(t3));
u3(1) = u0;
for i = 2:length(t3)
    u3(i) = u3(i-1) + h * f(u3(i-1));
    temp = inf;
    while abs(temp - u3(i)) > 1e-6
        temp = u3(i);
        u3(i) = u3(i-1) + h/2 * (f(u3(i-1)) + f(u3(i)));
    end
end
% 梯形法（方式三）
t4 = linspace(0, 1, 1/h);
u4 = zeros(size(t4));
u4(1) = u0;
for i = 2:length(t4)
    u4(i) = (1 + h/4)/(1 - h/4) * u4(i-1);
end
eps=2*1e-3;
for i=1:length(t4)
    if abs(u4(i)-u1(i))>eps
        break
    end
end
% 画图
figure;
% Euler法
subplot(3, 2, 1, 'Position', [0.06, 0.7, 0.4, 0.25]);
plot(t1, u1, 'k', rt, ru, 'r');
xlabel('X');
ylabel('Y');
title('Euler法');
legend('Euler法数值解', '真实解');
% 梯形法（方式一）
subplot(3, 2, 2, 'Position', [0.55, 0.7, 0.4, 0.25]);
plot(t2, u2, 'k', rt, ru, 'r');
xlabel('X');
ylabel('Y');
title('梯形法（方式一）');
legend('梯形法数值解', '真实解');
% 梯形法（方式二）
subplot(3, 2, 3, 'Position', [0.06, 0.35, 0.4, 0.24]);
plot(t3, u3, 'k', rt, ru, 'r');
xlabel('X');
ylabel('Y');
title('梯形法（方式二）');
legend('梯形法数值解', '真实解');
% 梯形法（方式三）
subplot(3, 2, 4, 'Position', [0.55, 0.35, 0.4, 0.24]);
plot(t4, u4, 'k', rt, ru, 'r');
xlabel('X');
ylabel('Y');
title('梯形法（方式三）');
legend('梯形法数值解', '真实解');

% 误差比较（对数尺度）
subplot(3, 2, 5, 'Position', [0.06, 0.07, 0.4, 0.18]);
plot(t1, abs(u1 - real_u(t1)), ...
     t2, abs(u2 - real_u(t2)), ...
     t3, abs(u3 - real_u(t3)), ...
     t4, abs(u4 - real_u(t4)));
set(gca, 'YScale', 'log');
xlabel('X');
ylabel('Y');
title('误差比较（对数尺度}');
legend('Euler法误差', '梯形法（方式一）误差', ...
       '梯形法（方式二）误差', '梯形法（方式三）误差');
% 步长与误差比较
subplot(3, 2, 6, 'Position', [0.55, 0.07, 0.4, 0.18]);
hold on;
plot(t1,abs(-u1+real_u(t1)),'b',t4,abs(-u4+real_u(t4)),'c',t1(i),abs(-u1(i)+real_u(t1(i))),'*',t4(i),abs(-u4(i)+real_u(t4(i))),'*')
xlabel('步长');
ylabel('误差');
title('误差相近时的步长与误差');
legend('Euler','梯形法')
set(gca, 'YScale', 'log');
set(gca,'XMinortick','on');
set(gca,'YMinortick','on');
hold off;
print('2', '-dpng');