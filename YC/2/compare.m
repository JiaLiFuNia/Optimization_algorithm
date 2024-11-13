clear;
close all;

name = '4_3';
n = 4;

load(['time',name,'.mat']);
load(['k',name,'.mat']);

set(gcf,'Position',[200,300,900,500])

subplot(1,2,1);
plot(1:n,log10(TIME),'LineWidth',1.2);
title('迭代时间 time 随规模的变化')
xlabel('规模n');
ylabel('$log_{10}(time)$','Interpreter','latex');
xticks(1:1:n) % x的范围
xticklabels({'10','100','1000','2000','5000'}) % 修改x的标签
legend({'基本牛顿法','牛顿-最速下降法','最速下降法','共轭梯度法','修正牛顿法',})

subplot(1,2,2);
plot(1:n,log10(K),'LineWidth',1.2);
ylim([0.5,4]);
title('迭代次数 k 随规模的变化')
xlabel('规模n');
ylabel('$log_{10}(k)$','Interpreter','latex');
xticks(1:1:n) % x的范围
xticklabels({'10','100','1000','2000','5000'}) % 修改x的标签
legend({'基本牛顿法','牛顿-最速下降法','最速下降法','共轭梯度法','修正牛顿法',})

print([name,'.png'],'-dpng','-r600');