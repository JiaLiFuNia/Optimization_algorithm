clear;
close all;
figure(1);

subplot(2,2,1);
x = -4*pi:0.01:4*pi;
f_sin = sin(x);
plot(x,f_sin,'b')
yticks(-1:0.2:1); % y的范围
xlim([x(1),x(end)]); % x的范围
xticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi]) % x的范围
xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'}) % 修改x的标签
yMinorTick = 'on'; % 次刻度
legend('Labels',{'f(x)=sin(x)'},'Orientation','horizontal'); % 图例
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
set(gca, 'Position',[0.05,0.55,0.4,0.4]);

subplot(2,2,2);
x = -2:0.1:2;
y = -2:0.1:2;
[X,Y] = meshgrid(x,y);
Z = sin(X).*Y;
[M,c] = contourf(Z);
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
clabel(M,c)
title('z=sin(x)*y 的等值线图');
set(gca, 'Position',[0.53,0.55,0.4,0.4]);

subplot(2,2,3);
mesh(X,Y,Z);
colormap winter;
colorbar;
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
zlabel('z','Interpreter','latex');
set(gca, 'Position',[0.04,0.05,0.4,0.4]);

subplot(2,2,4);
mesh(X,Y,Z);
colorbar;
hold on;
[M,c] = contourf(X,Y,Z);
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
zlabel('z','Interpreter','latex');
hold off;
set(gca, 'Position',[0.53,0.05,0.4,0.4]);

set(gcf,'Position',[200,300,1500,600]) % 设置画布大小
saveas(figure(1), 'example2', 'png');