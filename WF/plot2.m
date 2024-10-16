clear;
close;
figure(1);
x = -4*pi:0.001:4*pi;
f_sin = sin(x);
f_cos = cos(x);
f_3 = sin(x).*cos(x);
plot(x,f_sin,'b',x,f_cos,'r',x,f_3,'m');
xlim([-4*pi,4*pi])
legend('Labels',{'f(x)=sin(x)','g(x)=cos(x)','h(x)=sin(x)cos(x)'});
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
set(gca,'Position',[0.1 0.1 0.86 0.85],'YMinorTick','on')
gca.XAxis.MinorTickValues = -1:0.1:1;
% set(gcf,'Position',[200,300,1200,600]) % 设置画布大小
saveas(gcf, 'example1', 'png');