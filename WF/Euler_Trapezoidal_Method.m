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

%% figure(6)
[euler_u,tx_u] = figure6();

%% plot
close all;
set(gcf,'Position',[200,100,800,600])

subplot(3,2,1);
plot(t,euler_un);
hold on;
plot(t,true_un);
hold off;
xlim([0,1]);
title('Euler 法');
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xticks(0:0.1:1);
box on;
legend({'真解','数值解'});
set(gca,'Position',[0.06,0.72,0.43,0.25],'YMinorTick','on');

subplot(3,2,2);
plot(t,tx_un(:,1));
hold on;
plot(t,true_un);
hold off;
title('梯形法（方式1）');
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xticks(0:0.1:1);
box on;
legend({'真解','数值解'});
set(gca,'Position',[0.55,0.72,0.43,0.25],'YMinorTick','on');

subplot(3,2,3);
plot(t,tx_un(:,2));
hold on;
plot(t,true_un);
hold off;
title('梯形法（方式2）');
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xticks(0:0.1:1);
box on;
legend({'真解','数值解'});
set(gca,'Position',[0.06,0.39,0.43,0.24],'YMinorTick','on');

subplot(3,2,4);
plot(t,tx_un(:,3));
hold on;
plot(t,true_un);
hold off;
title('梯形法（方式3）');
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xticks(0:0.1:1);
box on;
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
box on;
legend({'梯形法1','梯形法2','梯形法3','Euler 法'});
set(gca,'Position',[0.06,0.06,0.43,0.24],'YMinorTick','on');

subplot(3,2,6);
for i = 1:size(euler_u,1)
    error = mean(abs(euler_u{i,1}-euler_u{i,2}));
    plot(euler_u{i,3},error,'b.');
    hold on;
end
for i = 1:size(tx_u,1)
    tx_ = tx_u{i,1};
    error = mean(abs(tx_(:,1)'-tx_u{i,2}));
    plot(tx_u{i,3},error,'r.');
    hold on;
end
hold off;
title('误差随步长的变化');
xlabel('$h$','Interpreter','latex');
ylabel('$error$','Interpreter','latex');
box on;
ax = gca;
ax.YAxis.Exponent = -3;
set(gca,'Position',[0.55,0.06,0.43,0.24],'YMinorTick','on');

saveas(figure(1), 'figure1', 'png');

%%
function [euler_u,tx_u] = figure6(~)
ut = @(t) exp(1/2*t);
ftu = @(t,u) 1/2*u;
u0 = 1;
Node = 1500;
% Euler
euler_u = {};
i = 0;
for node = 10:2:Node
    i = i + 1;
    t = linspace(0,1,node);
    h = 1/(node-1);
    euler_un(1) = u0;
    for k = 1:length(t)-1
        euler_un(k+1) = euler_un(k) + h*ftu(t(k),euler_un(k));
    end
    euler_u{i,1} = euler_un;
    euler_u{i,2} = ut(t);
    euler_u{i,3} = h;
end

% 梯形法
tx_u = {};
i = 0;
for node = 10:2:Node
    i = i + 1;
    t = linspace(0,1,node);
    h = 1/(node-1);
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
        while 1
            un_new = tx_un(k,2) + 0.5*h*(ftu(t(k),tx_un(k,2)) + ftu(t(k+1),tx_un(k+1,2)));
            if(abs(un_new - tx_un(k+1,2)) <= 1e-6)
                break;
            end
            tx_un(k+1,2) = un_new;
        end
    end

    % 3
    for k = 1:length(t) - 1
        tx_un(k+1,3) = tx_un(k,3)*(1+h/4)/(1-h/4);
    end

    tx_u{i,1} = tx_un;
    tx_u{i,2} = ut(t);
    tx_u{i,3} = h;
end

% 比较Euler法与梯形法显示公式的误差
% res = {};
% k = 1;
% for i = 1:size(tx_u,1)
%     tx_ = tx_u{i,1};
%     tx_error = mean(abs(tx_(:,3)'-tx_u{i,2}));
%     for m = 1:size(euler_u,1)
%         euler_error = mean(abs(euler_u{m,1}-euler_u{m,2}));
%         if abs(euler_error-tx_error) <= 1e-6
%             res{k,i} = [i,m];
%             k = k + 1;
%         end
%     end
% end
end