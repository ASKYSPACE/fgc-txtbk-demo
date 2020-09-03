clc
clear
close all
disp('%--------------------------------------------------------------%');
disp(' #=============================================================#');
disp(' |                         C-GO SCOP                           |');
disp(' | (Convex-programming based Guidance Optimization using SCOP) |');
disp(' |           (MATLAB R2016a + ecos  Full Edition)              |');
disp(' |                    Huan Jiang, Xinfu Liu                    |');
disp(' |     (Official site: http://www.asky.ac.cn/paper code)       |');
disp(' |  This is  COS LAB @ Beijing Institute of Tech.(03/21/2019)  |');
disp(' #==============================================================#');
g0 = 9.801;
emiu = 398600E9;      % Unit  m^3/s^2
Re  = 6378;           % Unit  m^3/s^2
hpf = 400.2;          % Unit  km
haf = 409.5;          % Unit  km
oba  = (haf+hpf)/2+Re;% Unit  km
obe  = (haf+Re)/oba-1;
flag = 1;    % 是否考虑轨道运动
Tall = 200; % 300s
N    = 101;
time = linspace(0,Tall,N);
tau  = time(3)-time(2);
flagobs  = 1 ; % 是否考虑避障
flagnote = 0;  % 是否记录结果数据

sate_pv01 = [-5;-35;0;-0.01;0.01;0];
sate_pv02 = [-10;-70;0;-0.01;0.01;0];
sate_pv03 = [0;-100;0;-0.02;0.01;0];
sate_pv0   = [sate_pv01 sate_pv02 sate_pv03];% postion/km velocity km/s
sate_R  = [6 6 6];
obaef  = [oba*1000,obe,0];%7378145;  4216400 e=0.3-0.7
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
sate_p = zeros(N,3,size(sate_pv0,2));
for i = 1:size(sate_pv0,2)
    sate_pv = obsevol(sate_pv0(:,i),time,obaef,emiu,opts,flag);
    sate_p(:,:,i) = sate_pv(:,1:3);
end
[dfenv,d2fenv,renv] = envpara(obaef,time,emiu);
%%
x0  = -10; % 0 -10 -20
y0  = -110; %-110
z0  = 0;
Vx0 = 0; 
Vy0 = 1;  
Vz0 = 0;
xf  = 0;
yf  = -10;
zf  = 0;
Vxf = 0;
Vyf = 0;  
Vzf = 0;
Icn = [0;-1;0];
alpha  = 10*pi/180;  %A half-angle 10-15 deg
m0 = 1000;%1905
Isp = 600/g0; % 发动机比冲 225s 喷射速率ve=g0*Isp  m/s
% Thrustmax = 16572.7;
Thrustmax = 15;
Speedmax   = 20;
Vmax = Speedmax;
damax = .0001;
roumax = 1.0*Thrustmax;
roumin = 0.2*Thrustmax;
% Csp = 1/(cos(27*pi/180)*g0*Isp);
Csp = 1/(g0*Isp);
lm0 = log(m0);
lmt0 = log(m0-Csp*roumax.*time);
lmt1 = log(m0-Csp*roumin.*time);
miu1 = roumin*exp(-lmt0);
miu2 = roumax*exp(-lmt0);
S  = 21;
%%
mk(1,:) = m0-Csp*roumax.*time;
Xk(1,:) = linspace(x0,xf,N);
Yk(1,:) = linspace(y0,yf,N);
Zk(1,:) = linspace(z0,zf,N);
vxk(1,:) = linspace(Vx0,Vxf,N);
vyk(1,:) = linspace(Vy0,Vyf,N);
vzk(1,:) = linspace(Vz0,Vzf,N);
uxk(1,:) = linspace(0,0,N-1);
uyk(1,:) = linspace(0,0,N-1);
uzk(1,:) = linspace(0,0,N-1);
tol = norm([x0-xf;y0-yf;z0-zf],2)*1.0e-4;% 0.01m
%%   %-----------------------------边界与约束----------------------%
for s = 1:S-1
    Xk_temp = Xk(s,:);Yk_temp = Yk(s,:);Zk_temp = Zk(s,:);
    for satei = 1:length(sate_R)
        for nn = 1:N
            Xt = [Xk_temp(nn);Yk_temp(nn);Zk_temp(nn)];
            p = sate_p(nn,:,satei)';
            Omega0(nn,satei) = sate_R(satei)^2 - (Xt - p)'*(Xt - p);
            DOmega0(:,nn,satei)  = -2*(Xt - p);
            lc = sate_R(satei)*(Xt-p)/norm(Xt-p);
            lmove(nn,satei) = ( Omega0(nn,satei) + DOmega0(:,nn,satei)'*(p+lc-Xt));
        end
    end
    
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
    F  = F +  [vx(1,1) == Vx0];
    F  = F +  [vx(1,N) == Vxf];
    F  = F +  [vy(1,1) == Vy0];
    F  = F +  [vy(1,N) == Vyf];
    F  = F +  [vz(1,1) == Vz0];
    F  = F +  [vz(1,N) == Vzf];
    %     F  = F +  [y <= yf];
    
    
    for i = 1:N-1
        % First-order Euler method
        F  = F + [x(1,i+1)==x(1,i)+tau*vx(1,i)];
        F  = F + [y(1,i+1)==y(1,i)+tau*vy(1,i)];
        F  = F + [z(1,i+1)==z(1,i)+tau*vz(1,i)];
        F  = F + [lm(1,i+1)==lm(1,i)-tau*Csp*sig(1,i)];
        F  = F + [vx(1,i+1)==vx(1,i)+tau*((dfenv(i)^2+2*emiu/renv(i)^3)*x(1,i)...
            +d2fenv(i)*y(1,i)+2*dfenv(i)*vy(1,i)+ux(1,i))];
        F  = F + [vy(1,i+1)==vy(1,i)+tau*((-d2fenv(i)^2)*x(1,i)...
            +(dfenv(i)^2-1*emiu/renv(i)^3)*y(1,i)+(-2*dfenv(i))*vx(1,i)+uy(1,i))];
        F  = F + [vz(1,i+1)==vz(1,i)+tau*((-1*emiu/renv(i)^3)*z(1,i)+uz(1,i))];
        
    end
    %         for i = 1:N-2
    %         F  = F + [norm([ux(1,i+1)-ux(1,i);uy(1,i+1)-uy(1,i);uz(1,i+1)-uz(1,i)],2)<=tau*damax ];
    %         end
    
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
    
    for i = 1:N
        F  = F + [norm([x(1,i);y(1,i);z(1,i)],2)<=...
            (x(1,i)*Icn(1)+y(1,i)*Icn(2)+z(1,i)*Icn(3))/cos(alpha)];
    end
    
    for j = 1:N-1
        F  = F + [norm([ux(1,j);uy(1,j);uz(1,j)],2)<=sig(1,j)];
        F  = F + [miu1(j)*(1-lm(1,j)+lmt0(j)+0.5*(lm(1,j)-lmt0(j))^2)<=sig(1,j)];
        F  = F + [sig(1,j)<=miu2(j)*(1-lm(1,j)+lmt0(j))];
        F  = F + [lmt0(j)<=lm(1,j)];
        F  = F + [lm(1,j)<=lmt1(j)];
    end
    
    if flagobs == 1
        for ocj = 1:length(sate_R)
            for oi = 1:N
                F  = [F; Omega0(oi,ocj)-lmove(oi,ocj)+DOmega0(:,oi,ocj)'...
                    *[(x(oi)-Xk_temp(oi));(y(oi)-Yk_temp(oi));(z(oi)-Zk_temp(oi))]<= 0];%
            end
        end
    end
    
    
    ops = sdpsettings('solver','mosek','verbose',1);%mosek ipopt ecos
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
    Fobj(s) = double(obj);
    mk(s+1,:) = m;
    Xk(s+1,:) = X; Yk(s+1,:) = Y; Zk(s+1,:) = Z;
    vxk(s+1,:) = Vx; vyk(s+1,:) = Vy; vzk(s+1,:) = Vz;
    uxk(s+1,:) = Ux; uyk(s+1,:) = Uy; uzk(s+1,:) = Uz;
    
%     lmt0 = log(m);
%     miu1 = roumin*exp(-lmt0);
%     miu2 = roumax*exp(-lmt0);
    
    errorx(s) = max(abs(Xk(s+1,:)-Xk(s,:)));
    errory(s) = max(abs(Yk(s+1,:)-Yk(s,:)));
    errorz(s) = max(abs(Zk(s+1,:)-Zk(s,:)));
    errorvx(s) = max(abs(vxk(s+1,:)-vxk(s,:)));
    errorvy(s) = max(abs(vyk(s+1,:)-vyk(s,:)));
    errorvz(s) = max(abs(vzk(s+1,:)-vzk(s,:)));    
    errorX(s,:)=[errorx(s),errory(s),errorz(s),errorvx(s),errorvy(s),errorvz(s)];
    fprintf('C-GO completed %d times optimization,the obj is %d.\n ',s,Fobj(s));
    fprintf('The iteration errors are: \n %d m, %d m, %d m,%d m/s, %d m/s, %d m/s.\n',errorX);
    if errorx(s) <= tol && errory(s) <= tol && errorz(s) <= tol &&...
            errorvx(s) <= tol && errorvy(s) <= tol && errorvz(s) <= tol
        break
    elseif flagobs == 0
        break
    end
    % >>>>>>> plot the iteration results <<<<<<<<< %
    lineW = 2;
    figure(1)
    hold on;
    subplot(2,1,1);
    plot((0:N-1)*tau,X,'-r', 'LineWidth', lineW); grid on; hold on;
    plot((0:N-1)*tau,Y,'-b', 'LineWidth', lineW); hold on;
    plot((0:N-1)*tau,Z,'-k', 'LineWidth', lineW); hold on;
    xlabel('time (s)'); ylabel('position (m)');
    subplot(2,1,2);
    hold on;
    plot((0:N-1)*tau,Vx,'-r', 'LineWidth', lineW); grid on; hold on;
    plot((0:N-1)*tau,Vy,'-b', 'LineWidth', lineW); hold on;
    plot((0:N-1)*tau,Vz,'-k', 'LineWidth', lineW); hold on;
    xlabel('time (s)'); ylabel('velocity (m/s)');
    figure(2)
    hold on;
    plot3(Y,X,Z,'b', 'LineWidth', lineW);axis equal;
    xlabel('postion y(m)');ylabel('position x(m)');zlabel('position z(m)');
    hold on;grid on;
    for i = 1:length(sate_R)
        plot3(sate_p(:,2,i),sate_p(:,1,i),sate_p(:,3,i),'r', 'LineWidth', lineW*2);
    end
    Vv = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
    figure(12)
    subplot(3,1,1)
    hold on;
    plot((0:N-1)*tau,Vx,'Linewidth',2);
    xlabel('time(s)');ylabel('Vx(m/s)');grid on;hold on;
    subplot(3,1,2)
    hold on;
    plot((0:N-1)*tau,Vy,'Linewidth',2);
    xlabel('time(s)');ylabel('Vy(m/s)');grid on;hold on;
    subplot(3,1,3)
    hold on;
    plot((0:N-1)*tau,Vz,'Linewidth',2);
    xlabel('time(s)');ylabel('Vz(m/s)');grid on;hold on;
    figure(13)
    hold on;
    plot((0:N-1)*tau,Vv,'Linewidth',2);
    xlabel('time(s)');ylabel('V(m/s)');grid on;hold on;
    figure(3);
    subplot(2,1,1);
    hold on;
    plot((0:N-2)*tau,sqrt(Tx.^2+Ty.^2+Tz.^2),'-r', 'LineWidth', lineW);
    xlabel('time (s)'); ylabel('throttle level'); grid on;
    subplot(2,1,2);
    hold on;
    plot((0:N-1)*tau,m, 'LineWidth', lineW);
    xlabel('time (s)'); ylabel('mass (kg)'); grid on;
    figure(5)
    hold on;
    plot(Y,X,'k', 'LineWidth', lineW);
    xlabel('y(m)');ylabel('x(m)');
    axis equal;
    set(gca,'XDir','reverse');
    %     axis image;
    hold on;grid on;
end

%% >>>>>>> plot the results <<<<<<<<< %%
clc
iters = size(Xk);
step = iters(1); 
date = datetime;
strtitle = strcat('RP',num2str(date.Month),num2str(date.Day),'_',...
    num2str(date.Hour),num2str(date.Minute),num2str(fix(date.Second)));
if flagnote == 1
    save(strtitle);
end
close all;
lineW = 2;
figure(1)
subplot(2,1,2);
hold on;
plot((0:N-1)*tau,X,'-r', 'LineWidth', lineW); grid on; hold on;
plot((0:N-1)*tau,Y,'--b', 'LineWidth', lineW); hold on;
plot((0:N-1)*tau,Z,'-.k', 'LineWidth', lineW); hold on;
xlabel('time (s)'); ylabel('position (m)');
legend('x','y','z');
subplot(2,1,1);
hold on;
plot3(Y,X,Z,'-k', 'LineWidth', lineW);
% axis equal;
set(gca,'XDir','reverse');
grid on;
xlabel('postion y(m)');ylabel('position x(m)');zlabel('position z(m)');
%%
figure(2)
hold on;
% axis([y0,yf,-40,30]);
plot3(Y,X,Z,'-k', 'LineWidth', lineW);
% axis equal;
xlabel('postion y(m)');ylabel('position x(m)');zlabel('position z(m)');
quiver3(Y(2:end),X(2:end),Z(2:end),Ty,Tx,Tz,1,'-k','Linewidth',0.5,'ShowArrowHead','off');
set(gca,'XDir','reverse');
hold on;grid on;
gamma = zeros(N-1,1);
for i = 1 :N-1
    U= [Tx(i),Ty(i),Tz(i)];
    V = [Vx(i+1),Vy(i+1),Vz(i+1)];
gamma(i) = acos(dot(U,V)/(norm(U)*norm(V)));
end
figure(2311)
subplot(211)
hold on;
plot((0:N-2)*tau,sqrt(Tx.^2+Ty.^2+Tz.^2),'-r', 'LineWidth', lineW);
xlabel('time (s)'); ylabel('throttle level (N)'); grid on;
subplot(212)
plot((0:N-2)*tau,Tx,'-r', 'LineWidth', lineW); grid on; hold on;
plot((0:N-2)*tau,Ty,'--b', 'LineWidth', lineW); hold on;
plot((0:N-2)*tau,Tz,'-.k', 'LineWidth', lineW); hold on;
xlabel('time (s)'); ylabel('throttle level (N)');
legend('Tx','Ty','Tz');
figure(15)
plot((0:N-2)*tau,gamma,':k','Linewidth',3);
xlabel('time (s)');ylabel('\gamma (deg)');grid on;hold on;
figure(21)
subplot(2,1,1)
hold on
plot((0:N-2)*tau,sqrt(Tx.^2+Ty.^2+Tz.^2),'-r', 'LineWidth', lineW);
xlabel('time (s)'); ylabel('throttle level (N)'); grid on;
subplot(2,1,2)
hold on;
plot((0:N-1)*tau,m,'k', 'LineWidth', lineW);
xlabel('time (s)'); ylabel('mass (kg)'); grid on;
%%
figure(211)
axis([y0,yf,-40,30]);
for i = 1:step
plot3(Yk(i,:),Xk(i,:),Zk(i,:),'--k', 'LineWidth', lineW*.5);
if i==step
 plot3(Yk(i,:),Xk(i,:),Zk(i,:),'k', 'LineWidth', lineW);     
end
% axis equal;
set(gca,'XDir','reverse');
hold on;grid on;
end
hold on;
% for i = 1:length(sate_R)
%  plot3(sate_p(:,2,i),sate_p(:,1,i),sate_p(:,3,i),'r', 'LineWidth', lineW);
%  hold on;
%  end
xlabel('postion y(m)');ylabel('position x(m)');zlabel('position z(m)');

%%
Vv = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
figure(212)
for i = 2:step
Vvk =sqrt(vxk(i,:).^2 + vyk(i,:).^2 + vzk(i,:).^2);
plot((0:N-1)*tau,Vvk,'--k', 'LineWidth', lineW*.5);
if i == step
plot((0:N-1)*tau,Vvk,'-k', 'LineWidth', lineW);    
end
xlabel('time (s)');ylabel('V (m/s)');grid on;hold on;
end
%%
figure(32)
hold on;
plot((0:N-1)*tau,Vv,'--k','Linewidth',2);
xlabel('time (s)');ylabel('V (m/s)');grid on;hold on;
figure(12)
subplot(2,1,2)
hold on;
plot((0:N-1)*tau,Vx,'-r','Linewidth',2);
xlabel('time (s)');ylabel('V (m/s)');grid on;hold on;
hold on;
plot((0:N-1)*tau,Vy,'--b','Linewidth',2);
xlabel('time(s)');ylabel('V (m/s)');grid on;hold on;
hold on;
plot((0:N-1)*tau,Vz,'-.k','Linewidth',2);
xlabel('time (s)');ylabel('V (m/s)');grid on;hold on;
legend('Vx','Vy','Vz');
subplot(2,1,1)
plot((0:N-1)*tau,Vv,'-k','Linewidth',2);
xlabel('time (s)');ylabel('V (m/s)');grid on;hold on;
%%
figure(3);
hold on;
plot((0:N-1)*tau,m,'k', 'LineWidth', lineW);
xlabel('time (s)'); ylabel('mass (kg)'); grid on;
figure(3131)
for i = 2:step
    Txk   = mk(i,2:end).*uxk(i,:);
    Tyk   = mk(i,2:end).*uyk(i,:);
    Tzk   = mk(i,2:end).*uzk(i,:);
    Tk =sqrt(Txk.^2 + Tyk.^2 + Tzk.^2);
plot((0:N-2)*tau,Tk,'--k', 'LineWidth', lineW*.5);
if i == step
plot((0:N-2)*tau,Tk,'k', 'LineWidth', lineW);
end
xlabel('time (s)');ylabel('throttle level (N)');grid on;hold on;
end
figure(312)
for i = 1:step
plot((0:N-1)*tau,mk(i,:),'--k', 'LineWidth', lineW*.5);
if i == step
plot((0:N-1)*tau,mk(i,:),'k', 'LineWidth', lineW);
end
xlabel('time (s)');ylabel('mass (kg)');grid on;hold on;
end
%%

figure(5)
distance = zeros(N,length(sate_R));
for j = 1:length(sate_R)
    distance(:,j) = sqrt((Y'-sate_p(:,2,j)).^2 + (X'-sate_p(:,1,j)).^2+ (Z'-sate_p(:,3,j)).^2)-sate_R(j);
    plot((0:N-1)*tau,distance(:,j),'k','LineWidth', lineW);
    hold on;grid on;
    xlabel('time (s)');ylabel('Distance (m)');
end
legend('一号障碍距离','二号障碍距离','三号障碍距离');
title('动态障碍与追踪航天器的距离');
%%

pic_num = 1;
theta=0:pi/100:2*pi;
figure(6)
hold on;
axis([y0,yf,-40,30]);
drawsate = fill([-1.4,-1,-1,-1.4],[0.2,0.2,0.45,0.45],'g');
satellite = fill([-1.4,-1,-1,-1.4],[0.2,0.2,0.45,0.45],'g');
drawtext=text(0,0,'');
gaptext = {'12','21'};
timetext=text(0,0,'');
for i = 1:N
    draw1=plot(Y(1:i),X(1:i),'k', 'LineWidth', lineW);
    xlabel('y(m)');ylabel('x(m)');grid on;
    %     axis equal;
    set(gca,'XDir','reverse');
    %     axis image;
    for j = 1:length(sate_R)
        drawsate(j)=fill(sate_p(1,2,j)+sate_R(j)*cos(theta),sate_p(1,1,j)+sate_R(j)*sin(theta),...
            [1 0.618 1]);%[0.618 0.618 0.618]
        gaptext = num2str(round(distance(i,j),2));
        gaptext = strcat('Distance=',gaptext,'m');
        drawtext(j)=text(sate_p(i,2,j)+6,sate_p(i,1,j)+6,gaptext);
        set(drawtext(j),'Position',[sate_p(i,2,j)+sate_R(j),sate_p(i,1,j)-1.2*sate_R(j)]);
        
        set(drawsate(j),'EdgeColor', 'r');
        set(drawsate(j),'XData',sate_p(i,2,j)+sate_R(j)*cos(theta));
        set(drawsate(j),'YData',sate_p(i,1,j)+sate_R(j)*sin(theta));
        satellite(j)=fill(sate_p(1,2,j)+1*cos(theta),sate_p(1,1,j)+1*sin(theta),'b');
        set(satellite(j),'XData',sate_p(i,2,j)+1*cos(theta));
        set(satellite(j),'YData',sate_p(i,1,j)+1*sin(theta));
        
        draw1=plot(sate_p(1:i,2,j),sate_p(1:i,1,j),'-k', 'LineWidth', 0.75*lineW);hold on;
        delete(timetext);
        timetext = title(strcat('Timer:',num2str(round((i-1)*tau,2)),'s'),...
            'HorizontalAlignment','center');%'FontSize',12,
    end
    pause(0.05);
    drawnow;
    F=getframe(gcf);
    if flagnote == 1
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,strcat(strtitle,'.gif'),'gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,strcat(strtitle,'.gif'),'WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;
    end
    if i ~= N
        delete(satellite);
        delete(drawsate);
        delete(drawtext);
    end
    hold on;
end
%%
figure(41)
semilogy(1:length(errorx),errorx,'-ok','LineWidth', lineW);
hold on;
semilogy(1:length(errory),errory,'--^b','LineWidth', lineW);
xlabel('iteration ');ylabel('error (m)');
grid on;title('位置迭代误差');
figure(42)
semilogy(1:length(errorvx),errorvx,'-ok','LineWidth', lineW);
hold on;
semilogy(1:length(errorvy),errorvy,'--^b','LineWidth', lineW);
xlabel('iteration ');ylabel('error (m/s)');
grid on;title('速度迭代误差');
fuelk = zeros(1,length(errorx)-1);
for i =1:length(errorx)
fuelk(i) = mk(i+1,1)-mk(i+1,end);
end
figure(44)
plot(1:length(fuelk),fuelk*1000,'--^k','LineWidth', lineW);
grid on;
% title('质量消耗迭代');
xlabel('iteration ');ylabel('m (g)');
figure(46)
semilogy(1:length(fuelk)-1,abs(diff(fuelk)),'--^b','LineWidth', lineW);
xlabel('iteration ');ylabel('m (kg)');
grid on;title('质量消耗迭代误差');