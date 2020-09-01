clc
clear
close all
g0 = 9.801;
emiu = 398600E9;    % Unit  m^3/s^2
R0  = 6378135;      % Unit  m
Lnorm = R0;
Vnorm = sqrt(R0*g0);
Tnorm = sqrt(R0/g0);
anorm = g0;
miunorm  = (Lnorm^3)/(Tnorm^2)/emiu;
flag = 1; % 是否考虑轨道运动
Tall = 1000; % 75s
N    = 100;
time = linspace(0,Tall,N)/Tnorm;
tau  = time(3)-time(2);
obs_s1 = [-10/Lnorm;-110/Lnorm;-10/Lnorm;0;0;0];
obs_s2 = [-200/Lnorm;-110/Lnorm;-50/Lnorm;0;0;0];
obs_s3 = [-300/Lnorm;-410/Lnorm;-260/Lnorm;0;0;0];
obs   = [obs_s1 obs_s2 obs_s3];% postion/km velocity km/s
obsr  = [100 100 100];
taraef  = [7378145/Lnorm,0.0,0];%7378145;  4216400 e=0.3-0.7 目标轨道参数
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
obss = zeros(N,3,size(obs,2));
for i = size(obs,2)
    obsp = obsevol(obs(:,i),time,taraef,1,opts,flag);
    obss(:,:,i) = obsp(:,1:3);
end
[dfenv,d2fenv,renv] = envpara(taraef,time,1);
%%
% x00  = 1500; y00  = 0; z00  = 2000;
% Vx00 = 0; Vy00 = 0;  Vz00 = 0;
% xff  = 0;yff  = 0; zff  = 0;
% Vxff = 0; Vyff = 0;  Vzff = 0;
x00  = -10; y00  = -110; z00  = 0;
Vx00 = 0; Vy00 = 0;  Vz00 = 0;
xff  = 0;yff  = -10; zff  = 0;
Vxff = 0; Vyff = 0;  Vzff = 0;
Icn = [0;-1;0];%/Lnorm
Ca  = 45*pi/180;
m0 = 1000;%1905
Isp = 60; % 发动机比冲 225s 喷射速率ve=g0*Isp  m/s
% Thrustmax = 16572.7;
Thrustmax = 15;
Speedmax   = 20;
Thrustnorm = m0*g0;
x0  = x00/Lnorm;xf  = xff/Lnorm;
y0  = y00/Lnorm;yf  = yff/Lnorm;
z0  = z00/Lnorm;zf  = zff/Lnorm;
Vx_0 = Vx00/Vnorm;Vx_f = Vxff/Vnorm;
Vy_0 = Vy00/Vnorm;Vy_f = Vxff/Vnorm;
Vz_0 = Vz00/Vnorm;Vz_f = Vxff/Vnorm;
Vmax = Speedmax/Vnorm;
amax = 9;
roumax = 0.8*Thrustmax/Thrustnorm;
roumin = 0.0*Thrustmax/Thrustnorm;
Csp = 1/(cos(27*pi/180)*g0*Isp/Vnorm);
lm0 = log(1);
lmt0 = log(1-Csp*roumax.*time);
lmt1 = log(1-Csp*roumin.*time);
miu1 = roumin*exp(-lmt0);
miu2 = roumax*exp(-lmt0);
S  = 2;
%%   %-----------------------------边界与约束----------------------%
for s = 1:S-1
    x  = sdpvar(1,N);
    y  = sdpvar(1,N);
    z  = sdpvar(1,N);
    lm = sdpvar(1,N);
    vx = sdpvar(1,N);
    vy = sdpvar(1,N);
    vz = sdpvar(1,N);
    sig= sdpvar(1,N-1);
    ux = sdpvar(1,N-1);
    uy = sdpvar(1,N-1);
    uz = sdpvar(1,N-1);
    F  = [ x(1,1) == x0];
    F  = F + [x(1,N) == xf];
    F  = F +  [y(1,1) == y0];
    F  = F +  [y(1,N) == yf];
    F  = F +  [z(1,1) == z0];
    F  = F +  [z(1,N) == zf];
    F  = F +  [lm(1,1) == lm0];
    F  = F +  [vx(1,1) == Vx_0];
    F  = F +  [vx(1,N) == Vx_f];
    F  = F +  [vy(1,1) == Vy_0];
    F  = F +  [vy(1,N) == Vy_f];
    F  = F +  [vz(1,1) == Vz_0];
    F  = F +  [vz(1,N) == Vz_f];
%     F  = F +  [y <= yf];    
    
    
    for i = 1:N-1
        % First-order Euler method
        F  = F + [x(1,i+1)==x(1,i)+tau*vx(1,i)];
        F  = F + [y(1,i+1)==y(1,i)+tau*vy(1,i)];
        F  = F + [z(1,i+1)==z(1,i)+tau*vz(1,i)];
        F  = F + [lm(1,i+1)==lm(1,i)-tau*Csp*sig(1,i)];
        F  = F + [vx(1,i+1)==vx(1,i)+tau*((dfenv(i)^2+2/renv(i)^3)*x(1,i)...
            +d2fenv(i)*y(1,i)+2*dfenv(i)*vy(1,i)+ux(1,i))];
        F  = F + [vy(1,i+1)==vy(1,i)+tau*((-d2fenv(i)^2)*x(1,i)...
            +(dfenv(i)^2-1/renv(i)^3)*y(1,i)+(-2*dfenv(i))*vx(1,i)+uy(1,i))];
        F  = F + [vz(1,i+1)==vz(1,i)+tau*((-1/renv(i)^3)*z(1,i)+uz(1,i))];
    end
    
    %         for i = 1:N-1
    %             % First-order Euler method
    %             F  = F + [x(1,i+1)==x(1,i)+tau*vx(1,i)];
    %             F  = F + [y(1,i+1)==y(1,i)+tau*vy(1,i)];
    %             F  = F + [z(1,i+1)==z(1,i)+tau*vz(1,i)];
    %             F  = F + [lm(1,i+1)==lm(1,i)-tau*Csp*sig(1,i)];
    %             F  = F + [vx(1,i+1)==vx(1,i)+tau*(ux(1,i))];
    %             F  = F + [vy(1,i+1)==vy(1,i)+tau*(uy(1,i))];
    %             F  = F + [vz(1,i+1)==vz(1,i)+tau*(uz(1,i))];
    %         end
    
%         for i = 1:N
%          F  = F + [norm([x(1,i);y(1,i);z(1,i)],2)<=...
%              (x(1,i)*Icn(1)+y(1,i)*Icn(2)+z(1,i)*Icn(3))/cos(Ca)];
%         end
    
    for j = 1:N-1
        F  = F + [norm([ux(1,j);uy(1,j);uz(1,j)],2)<=sig(1,j)];
        F  = F + [miu1(j)*(1-lm(1,j)+lmt0(j)+0.5*(lm(1,j)-lmt0(j))^2)<=sig(1,j)];
        F  = F + [sig(1,j)<=miu2(j)*(1-lm(1,j)+lmt0(j))];
        %         F  = F + [lm0t(j)<=lm(1,j)];
        %         F  = F + [lm(1,j)<=lm1t(j)];               
    end
    
    
    ops = sdpsettings('solver','ecos','verbose',1);%mosek ipopt ecos
    obj = sum(sig);
    solvesdp(F, obj, ops);
    X   = double(x);
    Y   = double(y);
    Z   = double(z);
    Lm   = double(lm);
    Vx   = double(vx);
    Vy   = double(vy);
    Vz   = double(vz);
    Ux   = double(ux);
    Uy   = double(uy);
    Uz   = double(uz);
    m    = exp(Lm);
    Tx   = m(1:end-1).*Ux;
    Ty   = m(1:end-1).*Uy;
    Tz   = m(1:end-1).*Uz;
    
    %% >>>>>>> plot the results <<<<<<<<< %%
    lineW = 2;
    figure(1)
    subplot(2,1,1);
    plot((0:N-1)*tau*Tnorm,X*Lnorm,'-r', 'LineWidth', lineW); grid on; hold on;
    plot((0:N-1)*tau*Tnorm,Y*Lnorm,'-b', 'LineWidth', lineW); hold on;
    plot((0:N-1)*tau*Tnorm,Z*Lnorm,'-k', 'LineWidth', lineW); hold on;
    xlabel('time (s)'); ylabel('position (m)');
    subplot(2,1,2);
    plot((0:N-1)*tau*Tnorm,Vx*Vnorm,'-r', 'LineWidth', lineW); grid on; hold on;
    plot((0:N-1)*tau*Tnorm,Vy*Vnorm,'-b', 'LineWidth', lineW); hold on;
    plot((0:N-1)*tau*Tnorm,Vz*Vnorm,'-k', 'LineWidth', lineW); hold on;
    xlabel('time (s)'); ylabel('velocity (m/s)');
    figure(2)
    plot3(X*Lnorm,Y*Lnorm,Z*Lnorm,'b', 'LineWidth', lineW);
    xlabel('x(m)');ylabel('y(m)');zlabel('z(m)');
    hold on;grid on;
    Vv = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
    figure(12)
    subplot(3,1,1)
    plot((0:N-1)*tau*Tnorm,Vx*Vnorm,'Linewidth',2);
    xlabel('time(s)');ylabel('Vx(m/s)');grid on;hold on;
    subplot(3,1,2)
    plot((0:N-1)*tau*Tnorm,Vy*Vnorm,'Linewidth',2);
    xlabel('time(s)');ylabel('Vy(m/s)');grid on;hold on;
    subplot(3,1,3)
    plot((0:N-1)*tau*Tnorm,Vz*Vnorm,'Linewidth',2);
    xlabel('time(s)');ylabel('Vz(m/s)');grid on;hold on;
    figure(13)
    plot((0:N-1)*tau*Tnorm,Vv*Vnorm,'Linewidth',2);
    xlabel('time(s)');ylabel('V(m/s)');grid on;hold on;
    figure(3);
    subplot(2,1,1);
    plot((0:N-2)*tau*Tnorm,sqrt(Tx.^2+Ty.^2+Tz.^2)*Thrustnorm,'-r', 'LineWidth', lineW);
    xlabel('time (s)'); ylabel('throttle level'); grid on;
    subplot(2,1,2);
    plot((0:N-1)*tau*Tnorm,m*m0, 'LineWidth', lineW);
    xlabel('time (s)'); ylabel('mass (kg)'); grid on;
    figure(5)
    plot(Y*Lnorm,X*Lnorm,'k', 'LineWidth', lineW);
    xlabel('y(m)');ylabel('x(m)');axis equal;
    hold on;grid on;
end
%%
% figure(1)
% plot3(obsp(:,1)*Lnorm,obsp(:,2)*Lnorm,obsp(:,3)*Lnorm);
% grid on;
% figure(2)
% plot(time*Tnorm,obsp(:,1)*Lnorm);