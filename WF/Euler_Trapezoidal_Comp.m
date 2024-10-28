function [euler_h,tx_h] = Euler_Trapezoidal_Comp(~)
ut = @(t) exp(1/2*t);
ftu = @(t,u) 1/2*u;
u0 = 1;
Node = 1000;

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
res = {};
k = 1;
for i = 1:size(tx_u,1)
    tx_ = tx_u{i,1};
    tx_error = mean(abs(tx_(:,2)'-tx_u{i,2}));
    for m = 1:size(euler_u,1)
        euler_error = mean(abs(euler_u{m,1}-euler_u{m,2}));
        if abs(euler_error-tx_error) <= 1e-6
            res{k} = [i,m];
            k = k + 1;
        end
    end
end

pos = res{11};
tx_h = tx_u{pos(1),3};
euler_h = euler_u{pos(2),3};
end