clear;
close all;
format long;

ut = @(t) sqrt(1+2*t);
ftu = @(t,u) u-2*t/u;
u0 = 1;

h = 0.1;
t = 0:h:1;
true_un = ut(t)';

%% 三阶
three_un = ones(length(t),1);

for k = 1:length(t)-1
    K1 = ftu(t(k),three_un(k));
    K2 = ftu(t(k)+1/2*h,three_un(k)+1/2*h*K1);
    K3 = ftu(t(k)+3/4*h,three_un(k)+3/4*h*K2);
    three_un(k+1) = three_un(k) + 1/9*h*(2*K1+3*K2+4*K3);
end

%% 四阶
four_un = ones(length(t),1);

for k = 1:length(t)-1
    K1 = ftu(t(k),four_un(k));
    K2 = ftu(t(k)+1/3*h,four_un(k)+1/3*h*K1);
    K3 = ftu(t(k)+2/3*h,four_un(k)-1/3*h*K1+h*K2);
    K4 = ftu(t(k)+h,four_un(k)+h*K1-h*K2+h*K3);
    four_un(k+1) = four_un(k) + 1/8*h*(K1+3*K2+3*K3+K4);
end

%% Adams q=3
adams_un = ones(length(t),1);
q = 3;
adams_un(1:q+1) = four_un(1:q+1);

for k = q+1:length(t)-1
    adams_un(k+1) = adams_un(k) + 1/24*h*(55*ftu(t(k),adams_un(k))-59*ftu(t(k-1),adams_un(k-1))+37*ftu(t(k-2),adams_un(k-2))-9*ftu(t(k-3),adams_un(k-3)));
end

%% PECE
pece_un = ones(length(t),1);
pece_un(1:4) = four_un(1:4);

for k = 1:length(t)-4
    pece_un(k+4) = pece_un(k+3) + 1/24*h*(55*ftu(t(k+3),pece_un(k+3))-59*ftu(t(k+2),pece_un(k+2))+37*ftu(t(k+1),pece_un(k+1))-9*ftu(t(k),pece_un(k)));
    fe = ftu(t(k+4),pece_un(k+4));
    while 1
        pece_un(k+4) = pece_un(k+3) + 1/24*h*(9*fe+19*ftu(t(k+3),pece_un(k+3))-5*ftu(t(k+2),pece_un(k+2))+ftu(t(k+1),pece_un(k+1)));
        fe_new = ftu(t(k+4),pece_un(k+4));
        if(abs(fe-fe_new) < 1e-6)
            break;
        end
        fe = fe_new;
    end
end

%% PMECME
% pmecme_un = ones(length(t),1);
% pmecme_un(1:4) = four_un(1:4);
% 
% for k = 1:length(t)-4
%     pmecme_un(k+4) = pmecme_un(k+3) + 1/24*h*(55*ftu(t(k+3),pmecme_un(k+3))-59*ftu(t(k+2),pmecme_un(k+2))+37*ftu(t(k+1),pmecme_un(k+1))-9*ftu(t(k),pmecme_un(k)));
%     um_3 = pmecme_un(k+3);
%     while 1
%         um = pmecme_un(k+4) + 251/270*(um_3-pmecme_un(k+3));
%         if(abs(um-pmecme_un(k+4)) < 1e-6)
%             break;
%         end
%     end
% end

%% figure
close all;
set(gcf,'Position',[200,100,800,600])

subplot(3,2,1);
yyaxis left
plot(t,three_un,'.-',t,true_un);
ylabel('$u$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
yyaxis right
plot(t,log10(abs(three_un-true_un)));
ylim([-12,-2]);
ylabel('$log_{10}(|u(t)-u_t|)$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
legend({'数值解','解析解','误差'},'Location','southeast');
title('三阶 Runge-Kutta 法');
set(gca,'Position',[0.06,0.72,0.38,0.24],'YMinorTick','on');

subplot(3,2,2);
yyaxis left
plot(t,four_un,'.-',t,true_un);
ylabel('$u$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
yyaxis right
plot(t,log10(abs(four_un-true_un)));
ylim([-12,-2]);
ylabel('$log_{10}(|u(t)-u_t|)$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
legend({'数值解','解析解','误差'},'Location','southeast');
title('四阶 Runge-Kutta 法');
set(gca,'Position',[0.56,0.72,0.38,0.24],'YMinorTick','on');

subplot(3,2,3);
yyaxis left
plot(t,adams_un,'.-',t,true_un);
ylabel('$u$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
yyaxis right
plot(t,log10(abs(adams_un-true_un)));
ylim([-12,-2]);
ylabel('$log_{10}(|u(t)-u_t|)$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
legend({'数值解','解析解','误差'},'Location','southeast');
title('Adams 显示格式');
set(gca,'Position',[0.06,0.39,0.38,0.24],'YMinorTick','on');

subplot(3,2,4);
yyaxis left
plot(t,pece_un,'.-',t,true_un);
ylabel('$u$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
yyaxis right
plot(t,log10(abs(pece_un-true_un)));
ylim([-12,-2]);
ylabel('$log_{10}(|u(t)-u_t|)$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
legend({'数值解','解析解','误差'},'Location','southeast');
title('PECE 格式');
set(gca,'Position',[0.56,0.39,0.38,0.24],'YMinorTick','on');

subplot(3,2,5);
yyaxis left
plot(t,three_un,'.-',t,true_un);
ylabel('$u$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
yyaxis right
plot(t,log10(abs(three_un-true_un)));
ylim([-12,-2]);
ylabel('$log_{10}(|u(t)-u_t|)$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
legend({'数值解','解析解','误差'},'Location','southeast');
title('PMCEMC 格式');
set(gca,'Position',[0.06,0.06,0.38,0.24],'YMinorTick','on');

subplot(3,2,6);
plot(t,log10(abs([three_un,four_un,adams_un,pece_un]-true_un)));
ylim([-12,-2]);
ylabel('$log_{10}(|u(t)-u_t|)$',"Interpreter","latex");
xlabel('$t$','Interpreter','latex');
legend({'三阶RK','四阶RK','Adams','PECE'},'Location','southeast','NumColumns',2);
set(gca,'Position',[0.56,0.06,0.38,0.24],'YMinorTick','on');

saveas(figure(1), 'figure2', 'png');