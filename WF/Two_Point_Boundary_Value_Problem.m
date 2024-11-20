clear;
format long;

f = @(x) sin(x) + cos(x);
ut = @(x) sin(x);

a = 0;b = pi;
alpha = 0;beta = 0;
p = 1;r = 1;q = 0;

ai = @(h,p,r) 2/h+1;
ci = @(h,p,r) 2/h-1;
bi = @(h,p,r) 4/h;
gi = @(h,x) 2*h*(sin(x)+cos(x));

%% h = 0.1pi
h_1 = pi*10^(-1);
x_1 = a:h_1:b;
N = length(x_1)-1;

d_a = -ai(h_1,p,r)*ones(N-1,1);
d_b = bi(h_1,p,r)*ones(N-1,1);
d_c = -ci(h_1,p,r)*ones(N-1,1);

% A = spdiags([d_a,d_b,d_c],-1:1,N-1,N-1);
y = gi(h_1,x_1(2:N)');

u_1 = [alpha;ZGF(d_a,d_b,d_c,y,N-1);beta];
true_u1 = ut(x_1);

%% h = pi*10^(-6)
h_2 = pi*10^(-6);
x_2 = a:h_2:b;
N = length(x_2)-1;

d_a = -ai(h_2,p,r)*ones(N-1,1);
d_b = bi(h_2,p,r)*ones(N-1,1);
d_c = -ci(h_2,p,r)*ones(N-1,1);

A = spdiags([d_a,d_b,d_c],-1:1,N-1,N-1);
y = gi(h_2,x_2(2:N)');

tic
u_2_1 = [alpha;A\y;beta];
t1 = toc;

tic
u_2_2 = [alpha;ZGF(d_a,d_b,d_c,y,N-1);beta];
t2 = toc;

true_u2 = ut(x_2);

%% plot
close all;

Label_H = 0.07;
Label_W = 0.07;
Weight = 0.42;
Height = 0.39;
P = [Label_W-0.01,2*Label_H+Height+0.03,Weight,Height;
    2*Label_W+Weight,2*Label_H+Height+0.03,Weight,Height;
    Label_W-0.01,Label_H,Weight,Height;
    2*Label_W+Weight,Label_H,Weight,Height];

set(gcf,'Position',[200,100,750,550]);

subplot(2,2,1)
plot(x_1,u_1,'s-','LineWidth',1.2);hold on;
fplot(ut,'LineWidth',1.2);hold off;
ylim([0,1.2]);
xlim([0,pi]);
xticks(0:0.2*pi:pi)
xticklabels({'0','0.2\pi','0.4\pi','0.6\pi','0.8\pi','\pi'})
legend({'数值解','解析解'});
xlabel('$x$','Interpreter','latex','FontSize',10);
ylabel('$u$','Interpreter','latex','FontSize',10);
title('$h=\pi\times10^{-1}$','Interpreter','latex');
set(gca,'position',[0.06 0.57 0.42 0.39],'YMinorTick','on','XMinorTick','on')

subplot(2,2,2)
plot(x_1(2:10),log10(abs(u_1(2:10)'-true_u1(2:10))),'LineWidth',1.2);
xlim([0,pi]);
xticks(0:0.2*pi:pi);
xticklabels({'0','0.2\pi','0.4\pi','0.6\pi','0.8\pi','\pi'})
xlabel('$x$','Interpreter','latex','FontSize',10);
ylabel('$log_{10}(|u-u_h|)$','Interpreter','latex','FontSize',10);
title('误差')
set(gca,'position',[0.56 0.57 0.42 0.39],'YMinorTick','on','XMinorTick','on')

subplot(2,2,3)
plot(x_2,u_2_1,'s-','LineWidth',1.2,'MarkerIndices',1:50000:N);hold on;
fplot(ut,'--','LineWidth',1.5);hold off;
ylim([0,1.2]);
xlim([0,pi]);
xticks(0:0.2*pi:pi)
xticklabels({'0','0.2\pi','0.4\pi','0.6\pi','0.8\pi','\pi'})
xlabel('$x$','Interpreter','latex','FontSize',10);
ylabel('$u$','Interpreter','latex','FontSize',10);
title('$h=\pi\times10^{-6}$','Interpreter','latex');
legend({'数值解','解析解'});
set(gca,'position',[0.06 0.07 0.42 0.39],'YMinorTick','on','XMinorTick','on')

subplot(2,2,4)
plot(x_2(2:end-1),log10(abs(u_2_1(2:end-1)'-true_u2(2:end-1))),'s-','LineWidth',1.2,'MarkerIndices',1:50000:N);hold on;
plot(x_2(2:end-1),log10(abs(u_2_2(2:end-1)'-true_u2(2:end-1))),'--','LineWidth',1.5);hold off;
xlim([0,pi]);
ylim([-16,-6]);
xticks(0:0.2*pi:pi)
xticklabels({'0','0.2\pi','0.4\pi','0.6\pi','0.8\pi','\pi'})
xlabel('$x$','Interpreter','latex','FontSize',10);
ylabel('$log_{10}(|u-u_h|)$','Interpreter','latex','FontSize',10);
title('误差');
legend({['求解器 t=',num2str(t1)],['追赶法 t=',num2str(t2)]});
set(gca,'position',[0.56 0.07 0.42 0.39],'YMinorTick','on','XMinorTick','on')

print('figure3.png','-dpng','-r600');
%% 追赶法
function x = ZGF(a,b,c,y,n)
x = zeros(n,1);
for i = 2:n
    w = a(i)/b(i-1);
    b(i) = b(i) - w*c(i-1);
    y(i) = y(i) - w*y(i-1);
end

x(n) = y(n)/b(n);
for i = n-1:-1:1
    x(i) = (y(i) - c(i)*x(i+1))/b(i);
end
end
