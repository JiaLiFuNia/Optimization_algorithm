close all;

name = '4_4';
n = 4;

load(['time',name,'.mat']);
load(['k',name,'.mat']);
N = [10,100,1000,2000];

X_pos = [0.08,0.56,0.08,0.56];
Y_pos = [0.57,0.57,0.08,0.08];

set(gcf,'Position',[200,300,900,500])
for i = 1:11:34
    subplot(2,2,mod(i,10));
    plot(log10(TIME(:,i:i+9)'))
    set(gca,'Position',[X_pos(mod(i,10)),Y_pos(mod(i,10)),0.4,0.36])
    ylim([-4,1]);
    title(['不同初值点下，运行时间的变化','(n=',num2str(N(mod(i,10))),')'])
    ylabel('$log_{10}(time)$','Interpreter','latex')
    xlabel('初值点')
    xticks(1:1:10) % x的范围
    xticklabels({'2.40','2.41','2.42','2.43','2.44','2.45','2.46','2.47','2.48','2.49'})
    legend({'基本牛顿法','牛顿-最速下降法','最速下降法','共轭梯度法','修正牛顿法'})
end
print([name,'.png'],'-dpng','-r600');

% set(gcf,'Position',[200,300,900,500])
% subplot(1,2,1);
% plot(1:n,log10(TIME),'LineWidth',1.2);
% title('迭代时间 time 随规模的变化')
% xlabel('规模n');
% ylabel('$log_{10}(time)$','Interpreter','latex');
% xticks(1:1:n) % x的范围
% xticklabels({'10','100','1000','2000','5000'}) % 修改x的标签
% legend({'基本牛顿法','牛顿-最速下降法','最速下降法','共轭梯度法','修正牛顿法',})
% 
% subplot(1,2,2);
% plot(1:n,log10(K),'LineWidth',1.2);
% ylim([0.5,4]);
% title('迭代次数 k 随规模的变化')
% xlabel('规模n');
% ylabel('$log_{10}(k)$','Interpreter','latex');
% xticks(1:1:n) % x的范围
% xticklabels({'10','100','1000','2000','5000'}) % 修改x的标签
% legend({'基本牛顿法','牛顿-最速下降法','最速下降法','共轭梯度法','修正牛顿法',})
% 
% print([name,'.png'],'-dpng','-r600');