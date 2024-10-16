clear;
close all;
format long;
FUN = {'@(x) exp(-x)+x.^2', '@(x) 3*x.^4-4*x.^3-12*x.^2', '@(x) x.^4+2*x+4', '@(x) x.^3-3*x+1'};

RESs = [];

for num = 1:length(FUN)
    phi = FUN{num};
    phi = str2func(phi);

    %% 进退法 十个初始点
    NUM = 0;
    RES = ones(8,10);
    for alpha0 = 0.5:0.1:1.4
        NUM = NUM + 1;
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
                    startLeftNode = min([Alpha,alpha(k+2)]);
                    startRightNode = max([Alpha,alpha(k+2)]);
                    break;
                end
            end
        end
        fprintf("单峰区间：[%f,%f]\n",startLeftNode,startRightNode);

        %% 0.618
        tic;
        a0 = startLeftNode;
        b0 = startRightNode;
        t = (sqrt(5)-1)/2;
        eps = 1e-5;
        delta = 0.15;

        a = [];
        b = [];
        p = [];
        q = [];

        a(1) = a0;
        b(1) = b0;
        p(1) = a(1) + (1-t)*(b(1)-a(1));
        q(1) = a(1) + t*(b(1)-a(1));
        h = abs(a(1) - b(1));

        i = 0;
        phip = phi(p(1));
        phiq = phi(q(1));
        fun_cal_count = 2;

        while(abs(q(i+1)-a(i+1)) > eps || abs(b(i+1)-p(i+1)) > eps || h > delta)
            if(phip <= phiq)
                a(i+2) = a(i+1);
                b(i+2) = q(i+1);
                phiq = phip;
                q(i+2) = p(i+1);
                p(i+2) = a(i+2) + (1-t)*(b(i+2)-a(i+2));
                phip = phi(p(i+2));
            else
                a(i+2) = p(i+1);
                b(i+2) = b(i+1);
                phip = phiq;
                p(i+2) = q(i+1);
                q(i+2) = a(i+2) + t*(b(i+2)-a(i+2));
                phiq = phi(q(i+2));
            end
            h = abs(a(i+2)-b(i+2));
            i = i + 1;
            fun_cal_count = fun_cal_count + 1;
        end
        if(phip <= phiq)
            minNode = p(i+1);
            minValue = phip;
        else
            minNode = q(i+1);
            minValue = phiq;
        end
        RES(1,NUM) = i;
        RES(3,NUM) = toc;
        RES(5,NUM) = fun_cal_count;
        RES(7,NUM) = minValue;
        disp('0.618法:');
        fprintf("极小点：%f, 最优解：%f, 迭代次数：%d\n",minNode,minValue,i);

        %% 抛物线法
        tic;
        h = abs(startLeftNode-startRightNode)/2;
        s0 = startLeftNode;
        s1 = startLeftNode + h;
        s2 = startRightNode;
        phi0 = phi(s0);
        phi1 = phi(s1);
        phi2 = phi(s2);

        eps = 1e-5;
        count = 0;
        fun_cal_count = 3;
        while(abs(s1-s0) >= eps && abs(s2-s1) >= eps)
            q_up = (s1^2-s2^2)*phi0+(s2^2-s0^2)*phi1+(s0^2-s1^2)*phi2;
            q_down = (s1-s2)*phi0+(s2-s0)*phi1+(s0-s1)*phi2;
            s_ba = 0.5 * q_up / q_down;
            phi_ba = phi(s_ba);
            fun_cal_count = fun_cal_count + 1;

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
        RES(2,NUM) = count;
        RES(4,NUM) = toc;
        RES(6,NUM) = fun_cal_count;
        RES(8,NUM) = phi(s1);
        disp('抛物线法：')
        fprintf("极小点：%f, 最优解：%f, 迭代次数：%d\n\n", s1,phi(s1),count);
    
    end
    RESs = [RESs;zeros(1,10);RES];
end
disp(RESs)