%% 进退法
clear;
close all;
format long;
phi = @(x) x.^2-6.*x+2;
subplot(1,2,1);
fplot(phi);
alpha0 = 1;
h0 = 0.1;
k = 0;

h = [];
phik = [];
alpha = [];

alpha(1) = alpha0;
phik(1) = phi(alpha0);
h(1) = h0;
Alpha = 0;

while(1)
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

%% 抛物线法
h = (b-a)/2;
s0 = a;
s1 = a + h;
s2 = b;
phi0 = phi(s0);
phi1 = phi(s1);
phi2 = phi(s2);

eps = 1e-5;
count = 0;
while(abs(s2-s0) >= eps && abs(s2-s1) >= eps)
    q_up = (s1^2-s2^2)*phi0+(s2^2-s0^2)*phi1+(s0^2-s1^2)*phi2;
    q_down = (s1-s2)*phi0+(s2-s0)*phi1+(s0-s1)*phi2;
    s_ba = 0.5 * q_up / q_down;
    phi_ba = phi(s_ba);
    if(phi1 <= phi_ba)
        if(s1 < s_ba)
            s2 = s_ba;
            phi2 = phi_ba;
        else
            s0 = s_ba;
            phi0 = phi_ba;
        end
    else
        if(s1 > s_ba)
            s2 = s1;
            s1 = s_ba;
            phi2 = phi1;
            phi1 = phi_ba;
        else
            s0 = s1;
            s1 = s_ba;
            phi0 = phi1;
            phi1 = phi_ba;
        end
    end
    count = count + 1;
end
fprintf("极小点：%f，最优解：%f, 迭代次数：%d\n", s1,phi(s1),count);
