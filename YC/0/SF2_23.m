clc;
clear;
phi = @(x) x.^2-6.*x+2;
alpha0 = 1;
h0 = 0.1;

alpha = alpha0;
alphak = alpha0;
hk = h0;
alphak1 = alphak + hk;

phik = feval(phi,alphak);
phik1 = feval(phi,alphak1);
k=0;
while((phik1<phik)||(k==0))
    if(phik1<phik)
        hk = 2*hk;
        alpha = alphak;
        alphak = alphak1;
        phik = phik1;
        k = k + 1;

        alphak1 = alphak + hk;
        phik1 = feval(phi, alphak1);
    else if(k==0)
            hk = - hk;
            alpha = alphak1;
            alphak = alphak;
            phik = phik;
            k=1;

            alphak1 = alphak + hk;
            phik1 = feval(phi, alphak1);
    end
    end
end
s0 = min(alpha,alphak1);
s2 = max(alpha,alphak1);
s1 = alphak;
eps = 1e-5;

phis0 = feval(phi,s0);
phis1 = feval(phi,s1);
phis2 = feval(phi,s2);
k = 0;
while (abs(s1 - s0) >= eps) && (abs(s2 - s1) >= eps)
    part1 = (s1^2-s2^2)*phis0+(s2^2-s0^2)*phis1+(s0^2-s1^2)*phis2;
    part2 = (s1-s2)*phis0+(s2-s0)*phis1+(s0-s1)*phis2;
    s_ba = 0.5 * part1 / part2;
    phis_ba = feval(phi,s_ba);
    if (phis1 <= phis_ba)
        if (s1 < s_ba)
            s2 = s_ba;
            phis2 = phis_ba;
        else
            s0 = s_ba;
            phis0 = phis_ba;
        end
    else
        if (s1 > s_ba)
            s2 = s1;
            s1 = s_ba;
            phis2 = phis1;
            phis1 = phis_ba;
        else
            s0 = s1;
            s1 = s_ba;
            phis0 = phis1;
            phis1 = phis_ba;
        end
    end
    k=k+1;
end
s1