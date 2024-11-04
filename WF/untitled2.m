clc,clear;%涕型法1
u(1)=1;x0=0;
h=0.1;
for i=1:10
    f=u(i)/2;
    u(i+1)=u(i)+f*h;
    for k=1:3
        u(i+1)=u(i)+h*(f+u(i+1)/2)/2;
    end
end
%%真姐
fun=@(t)exp(t/2);
t=0:0.1:1;
g=fun(t);
d(1,:)=u;
d(2,:)=g;
d(3,:)=t;
subplot(2,3,2)
plot(t,u,t,g,'.')
xlabel('x')
ylabel('y')
legend('数值解','真解')
title('梯形法1')
grid on
%%
clc,clear;%涕型法2
u(1)=1;x0=0;
h=0.1;
for i=1:10
    f=u(i)/2;
    u(i+1)=u(i)+f*h;
    %for k=1:3
    u(i+1)=u(i)+h*(f+u(i+1)/2)/2;%矫正
    %end
end
%%真姐
fun=@(t)exp(t/2);
t=0:0.1:1;
g=fun(t);
d(1,:)=u;
d(2,:)=g;
d(3,:)=t;
subplot(2,3,3)
plot(t,u,t,g,'.')
xlabel('x')
ylabel('y')
legend('数值解','真洁')
title('梯形法2')
grid on
%%
clc,clear;%涕型法3
u(1)=1;x0=0;
h=0.1;
for i=1:10
    f=u(i)/2;
    u(i+1)=u(i)+f*h;
    u(i+1)=(1+h/4)*u(i)/(1-h/4);
end
%%真姐
fun=@(t)exp(t/2);
t=0:0.1:1;
g=fun(t);
d(1,:)=u;
d(2,:)=g;
d(3,:)=t;
subplot(2,3,4)
plot(t,u,t,g,'.')
xlabel('x')
ylabel('y')
legend('数值解','真解')
title('梯形法3')
grid on
%%
%%Euler法'
clc,clear;
u(1)=1;
x0=0;
h=0.1;
for i=1:10
    f=u(i)/2;
    u(i+1)=u(i)+f*h;
end
%%真姐
fun=@(t)exp(t/2);
t=0:0.1:1;
g=fun(t);
subplot(2,3,1)
plot(t,u,t,g,'.')
xlabel('x')
ylabel('y')
legend('数值解','真解')
title('Euler法')
grid on
%%
clc,clear;
u1(1)=1;x0=0;
h=0.1;
for i=1:10
    f=u1(i)/2;
    u1(i+1)=u1(i)+f*h;
end
%%真姐
fun=@(t)exp(t/2);
t=0:0.1:1;
g=fun(t);
for i=1:11
    s1(i)=abs(log(u1(i)-g(i)));
end
%%
u2(1)=1;x0=0;
h=0.1;
for i=1:10
    f=u2(i)/2;
    u2(i+1)=u2(i)+f*h;
    for k=1:3
        u2(i+1)=u2(i)+h*(f+u2(i+1)/2)/2;
    end
end
%%真姐
fun=@(t)exp(t/2);
t=0:0.1:1;
g=fun(t);
for i=1:11
    s2(i)=abs(log(u2(i)-g(i)));
end
%%
u3(1)=1;x0=0;
h=0.1;
for i=1:10
    f=u3(i)/2;
    u3(i+1)=u3(i)+f*h;
    %for k=1:3
    u3(i+1)=u3(i)+h*(f+u3(i+1)/2)/2;
    %end
end
%%真姐
fun=@(t)exp(t/2);
t=0:0.1:1;
g=fun(t);
for i=1:11
    s3(i)=abs(log(u3(i)-g(i)));
end
%%
u4(1)=1;x0=0;
h=0.1;
for i=1:10
    f=u4(i)/2;
    u4(i+1)=u4(i)+f*h;
    u4(i+1)=(1+h/4)*u4(i)/(1-h/4);
    %end
end
%%真姐
fun=@(t)exp(t/2);
t=0:0.1:1;
g=fun(t);
for i=1:11
    s4(i)=abs(log(u4(i)-g(i)));
end
subplot(2,3,5)
plot(t,s1,t,s2,t,s3,t,s4);
legend('欧拉法','1','2','3')