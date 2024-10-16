% 进退法
clear;
close all;
phi = @(x) x.^2-6*x+2;
subplot(1,2,1);
fplot(phi);
alpha0 = 1;
h0 = 2;
a = 0; b = 1;
xlim([a,b]);
k = 0;

h = [];
phik = [];
alpha = [];

alpha(1) = alpha0;
phik(1) = phi(alpha0);
h(1) = h0;
Alpha = 0;


while 1
    alpha(k+2) = alpha(k+1) + h(k+1);
    phik(k+2) = phi(alpha(k+2));
    if(phik(k+2) < phik(k+1))
        h(k+2) = 2*h(k+1);
        Alpha = alpha(k+1);
        alpha(k+1) = alpha(k+2);
        phik(k+1) = phik(k+2);
        k = k + 1;
    else
        if(k == 0)
            h(2) = -h(1);
            Alpha = alpha(2);
            alpha(2) = alpha(1);
            phik(2) = phik(1);
            k = 1;
        else
            a = min([Alpha,alpha(k+2)]);
            b = max([Alpha,alpha(k+2)]);
            break;
        end
    end
end
subplot(1,2,2);
fplot(phi);
xlim([a,b]);
