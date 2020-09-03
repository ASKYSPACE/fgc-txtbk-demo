%Xm,Ym,Xt,Yt        飞行器的初始位置与终点位置
%theta0，thetaf     飞行器的初始以及终点轨迹角
%V0                 飞行器飞行速度
%tfmin，tfmax       飞行器最大飞行时长
%m                  飞行器质量
%S                  飞行器参考面积
%rou                空气密度
%g                  重力常数
%Cd                 阻力系数
%h                  迭代步长

%%
%设定初始位置，末端位置，飞行时间，末端路径角
clc;clear;
Xm =    0;Ym = 2500;
Xt = 5000;Yt = 0;
theta0 = 0; V0 = 300; 
thetaf = -45/57.3;
 tfmin = 25;t0 = 0;
 tfmax = 30;
 %%
 %设定参考面积，空气密度，质量，阻力系数，离散时间
global m S Cd g rou  
m=200;
S = 0.04;
rou = 1;
g=9.8;
Cd = 0.1;
%%
%给出一条初始零弹道
h = 0.2;
N = tfmax/0.2;
tic
U(1,:)=zeros(1,N);  %初始解猜想
[t,Data] = jifen(@Dynamics,[V0;theta0;Xm;Ym],t0,tfmax,N ,U(1,:));
Fu(:,1) = [Data(2,end)-thetaf;Data(3,end)-Xt;Data(4,end)-Yt;];
range = sqrt(Fu(2,1)^2+Fu(3,1)^2);
k = 1; 
J = 100;
%%
%开始迭代求解参考书中求解公式
while abs(J)>0.1
    for i = 1:N
        pGu = [0;h/Data(1,i);0;0;];
        for j = i:N
            pGx = [1-rou*S*Cd/m*Data(1,j)*h -g*cos(Data(2,j))*h 0 0;
                   (g*cos(Data(2,j))-U(k,j))/(Data(1,j)^2)*h 1+g*sin(Data(2,j))/Data(1,j)*h 0 0;
                   cos(Data(2,j))*h -(Data(1,j))*sin(Data(2,j))*h 1 0;
                   sin(Data(2,j))*h (Data(1,j))*cos(Data(2,j))*h 0 1;];
            pGu = pGx*pGu;
        end
        pYx = [0 1 0 0;0 0 1 0;0 0 0 1];
        pFu = pYx*pGu;
        dFu(:,i)=pFu;
    end
    s(:,k) = -(dFu')*pinv(dFu*(dFu'))*(dFu*(U(k,:)')-Fu(:,k))+U(k,:)';
    U(k+1,:) = U(k,:)-s(:,k)';
    k = k+1;
    [t,Data] = jifen(@Dynamics,[V0;theta0;Xm;Ym],t0,tfmax,N ,U(k,:));
    Fu(:,k) = [Data(2,end)-thetaf;Data(3,end)-Xt;Data(4,end)-Yt;];
    range = sqrt(Fu(2,k)^2+Fu(3,k)^2);
    J = norm(U(k,:))-norm(U(k-1,:));
    Jmin = 0.5*norm(U(k,:));
    if k>100
        break;
    end
    
end
toc
%%
%绘图
figure(1);
for i=2:k
    plot(t(1:N),U(i,:),'LineWidth',1.5);hold on;grid on;
end
legend('第1次迭代','第2次迭代','第3次迭代','第4次迭代','第5次迭代','第6次迭代','第7次迭代');
xlabel('T/s');ylabel('控制量a');
figure(2);
plot(t(1:N),U(k,:),'k','LineWidth',1.5);hold on;grid on;
xlabel('T/s');ylabel('控制量a');