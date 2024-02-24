% Zane Grothe
% AERO 7970
% HW 1
% 8/26/22

clear all
close all
clc

%% Problem 1 ~~~~~~~~~~~~~~~~~~~~

% Part b

b=[-10^2,-1,0,1,10^2,10^4];
x=linspace(0,1);

c=1;
while c < 7
    a=b(1,c);
    if a==1
        disp(sprintf('In Problem 1 Part b there is no solution for a=%.0f',a));
    else
        figure(c)
        P=exp(x*(a-1))/(a-1);
        Q=1-1/(a-1);
        R=exp(a*x);
        u=P+Q./R;
        plot(x,u)
        title(sprintf('Problem 1 Part b Graph for a=%.0f',a))
    end
    c=c+1;
end


%% Problem 2 ~~~~~~~~~~~~~~~~~~~~

% Part b

K=0.05;
A=1;
B=1;
t=linspace(0,100);

y=(1./exp(K*t)).*(exp(K*t).*(A*(144*K^2+pi^2)+144*B*K^2*sin(pi*t./12)...
    -12*pi*B*K*cos(pi*t./12)))./(K*(144*K^2+pi^2))...
    -((A*(144*K^2+pi^2)-12*pi*B*K)/(K*(144*K^2+pi^2)));

figure(7)
plot(t,y)
title('Blood Hormone Levels (BHL) vs Time')
xlabel('Time')
ylabel('BHL')
