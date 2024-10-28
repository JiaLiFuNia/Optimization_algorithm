clear;
close all;
format long;

ut = @(t) exp(1/2*t);
ftu = @(t,u) 1/2*u;
u0 = 1;

h = 0.1;
t = 0:h:1;
true_un = ut(t);

%% Euler
euler_un = ones(length(t),1);
euler_un(1) = u0;
for k = 1:length(t)-1
    euler_un(k+1) = euler_un(k) + h*ftu(t(k),euler_un(k));
end

%% 梯形法
tx_un = [];
tx_un(1,1:3) = u0;

% 1
for k = 1:length(t) - 1
    tx_un(k+1,1) = tx_un(k,1) + h*ftu(t(k),tx_un(k,1));
    tx_un(k+1,1) = tx_un(k,1) + 0.5*h*(ftu(t(k),tx_un(k,1)) + ftu(t(k+1),tx_un(k+1,1)));
end

% 2
for k = 1:length(t) - 1
    tx_un(k+1,2) = tx_un(k,2) + h*ftu(t(k),tx_un(k,2));
    i = 0;
    while 1
        un_new = tx_un(k,2) + 0.5*h*(ftu(t(k),tx_un(k,2)) + ftu(t(k+1),tx_un(k+1,2)));
        if(abs(un_new - tx_un(k+1,2)) <= 1e-6)
            break;
        end
        tx_un(k+1,2) = un_new;
        i = i + 1;
    end
end

% 3
for k = 1:length(t) - 1
    tx_un(k+1,3) = (1+h/4)/(1-h/4)*tx_un(k,3);
end

%% plot
close all;
set(gcf,'Position',[200,100,800,600])

subplot(3,2,1);
plot(t,euler_un);
hold on;
plot(t,true_un,'.-');
hold off;
xlim([0,1]);
title('Euler 法');
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xticks(0:0.1:1);
legend({'真解','数值解'});
set(gca,'Position',[0.06,0.72,0.43,0.25],'YMinorTick','on');

subplot(3,2,2);
plot(t,tx_un(:,1));
hold on;
plot(t,true_un,'.-');
hold off;
title('梯形法（方式1）');
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xticks(0:0.1:1);
legend({'真解','数值解'});
set(gca,'Position',[0.55,0.72,0.43,0.25],'YMinorTick','on');

subplot(3,2,3);
plot(t,tx_un(:,2));
hold on;
plot(t,true_un,'.-');
hold off;
title('梯形法（方式2）');
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xticks(0:0.1:1);
legend({'真解','数值解'});
set(gca,'Position',[0.06,0.39,0.43,0.24],'YMinorTick','on');

subplot(3,2,4);
plot(t,tx_un(:,3));
hold on;
plot(t,true_un,'.-');
hold off;
title('梯形法（方式3）');
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xticks(0:0.1:1);
legend({'真解','数值解'});
set(gca,'Position',[0.55,0.39,0.43,0.24],'YMinorTick','on');

subplot(3,2,5);
plot(t,log10(abs(tx_un-true_un')));
hold on;
plot(t,log10(abs(euler_un-true_un')));
hold off;
title('Euler 法与梯形法误差曲线');
xlabel('$t$','Interpreter','latex');
ylabel('$log_{10}(|u(t)-u_t|)$','Interpreter','latex');
xticks(0:0.1:1);
legend({'梯形法1','梯形法2','梯形法3','Euler 法'});
set(gca,'Position',[0.06,0.06,0.43,0.24],'YMinorTick','on');

%% figure(6)
% 当Euler法与梯形法误差很小时
[euler_h,tx_h] = Euler_Trapezoidal_Comp()

% Euler
h = euler_h;
euler_t = 0:h:1;
euler_true_un = ut(euler_t);

euler_un = ones(length(euler_t),1);
euler_un(1) = u0;
for k = 1:length(euler_t)-1
    euler_un(k+1) = euler_un(k) + h*ftu(euler_t(k),euler_un(k));
end

% 梯形法
h = tx_h;
t = 0:h:1;
tx_true_un = ut(t);

tx_un = [];
tx_un(1) = u0;

% 2
for k = 1:length(t) - 1
    tx_un(k+1) = tx_un(k) + h*ftu(t(k),tx_un(k));
    i = 0;
    while 1
        un_new = tx_un(k) + 0.5*h*(ftu(t(k),tx_un(k)) + ftu(t(k+1),tx_un(k+1)));
        if(abs(un_new - tx_un(k+1)) <= 1e-6)
            break;
        end
        tx_un(k+1) = un_new;
        i = i + 1;
    end
end

%% plot6
subplot(3,2,6);
plot(t,log10(abs(tx_un-tx_true_un)),'LineWidth',1.2);
hold on;
plot(euler_t,log10(abs(euler_un-euler_true_un')));
hold off;
title('Euler 法与梯形法的误差曲线');
xlabel('$t$','Interpreter','latex');
ylabel('$log_{10}(|u(t)-u_t|)$','Interpreter','latex');
xticks(0:0.1:1);
legend({['梯形法 (h=',num2str(tx_h),')'],['Euler 法 (h=',num2str(euler_h),')']});
set(gca,'Position',[0.55,0.06,0.43,0.24],'YMinorTick','on');

%% save
saveas(figure(1), 'figure1', 'png');