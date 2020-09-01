clc;
clear all;
close all;
format long;
disp('%--------------------------------------------------------------%');
disp(' #=============================================================#');
disp(' |                         C-GO SCOP                           |');
disp(' | (Convex-programming based Guidance Optimization using SCOP) |');
disp(' |           (MATLAB R2016a + ecos  Full Edition)              |');
disp(' |                    Huan Jiang, Xinfu Liu                    |');
disp(' |     (Official site: http://www.asky.ac.cn/paper code)       |');
disp(' |  This is  COS LAB @ Beijing Institute of Tech.(03/21/2019)  |');
disp(' #==============================================================#');


%% %------------------Modeling and Constraint setting------------%
disp('Modeling and Constraint Setting ...');
x00  = 0; y00  = 0;  z00  = 0;
xff  = 60;yff  = 60; zff  = 0;
Vx00 = 5; Vy00 = 3;  Vz00 = 0;
Vxff = 5; Vyff = 3;  Vzff = 0;
m0   = 1000;
z0   =log(m0);
Isp  = 100;
%====================[-SPECS-DJI-MAVIC2 S-Model]=========================%
rho_max = 20;            % Max Speed (near sea level)
rho_min = 10;            % Max level Speed 
Des_max   = 3;             % Max Descent Speed
As_max    = 5;             % Max Ascent Speed
Accel_max = 12;            % Max Acceleration
pullz_min = .05*rho_max; % Approximate minimum acceleration
Tilt_max  = 35;            % Max Tilt Angle
%====================[-SPECS-DJI-MAVIC2 S-Model]=========================%

flagobs = 0;
% method = 'Linearization';
method = 'Projection';
N   = 100;  % Discretization points
tau = 1/N;
g = 9.82;
Lnorm = norm([xff-x00;yff-y00]);
Vnorm = sqrt(g*Lnorm);
anorm = g;
Tnorm = Lnorm/Vnorm;

x0  = x00/Lnorm;
xf  = xff/Lnorm;
y0  = y00/Lnorm;
yf  = yff/Lnorm;
z0  = z00/Lnorm;
zf  = zff/Lnorm;
Vx_0 = V0*cos(psi0)*cos(phi0)/Vnorm;
Vy_0 = V0*cos(psi0)*sin(phi0)/Vnorm;
Vz_0 = V0*sin(psi0)/Vnorm;
Vx_f = Vf*cos(psif)*cos(phif)/Vnorm;
Vy_f = Vf*cos(psif)*sin(phif)/Vnorm;
Vz_f = Vf*sin(psif)/Vnorm;
Speedmax = Speed_max/Vnorm;
Vlmax   = Vl_max/Vnorm;
Vzmin   = -Des_max/Vnorm;
Vzmax   = As_max/Vnorm;
amax    = Accel_max/anorm;%12
uzmin   = pullz_min/anorm;
Tiltmax = Tilt_max*pi/180;
s=1;
tol = 1.0e-4;
Roc = R_oc;
obsp = obsevol(obs_p0,obs_v0,aef,time,flag);
[sit,dsit,r]=orbitevol(aef,time);
%%
S = 10;
S=S+1;
Lambda(1,:) = 0;
T1k(1) = sqrt((xf-x0)^2+(yf-y0)^2)/(0.5*(V0+Vf))/(N*tau);%
Xk(1,:) = linspace(x0,xf,N);
Yk(1,:) = linspace(y0,yf,N);
Zk(1,:) = linspace(z0,zf,N);
vxk(1,:) = linspace(Vx_0,Vx_f,N-1);
vyk(1,:) = linspace(Vy_0,Vy_f,N-1);
vzk(1,:) = linspace(Vz_0,Vz_f,N-1);
uxk(1,:) = linspace(1,1,N-2)*0.001;
uyk(1,:) = linspace(1,1,N-2)*0.001;
uzk(1,:) = linspace(1,1,N-2)*0.001;
F = [];
%%   %-----------------------------±ß½çÓëÔ¼Êø----------------------%
for s = 1:S-1
    Xk_temp = Xk(s,:);Yk_temp = Yk(s,:);Zk_temp = Zk(s,:);
    for oci = 1:length(Roc)
        for nn = 1:N
            for ti = 1:N
            Omega0(nn,oci)  = [(Xk_temp(nn) - obsp(oci));(Yk_temp(nn) - Yoc(oci))]'*inv(Aoc{oci,1})'/(Aoc{oci,1})*[(Xk_temp(nn) - Xoc(oci));(Yk_temp(nn) - Yoc(oci))]  - Roc(oci)^2;
            Dgc(:,nn,oci)  = 2*inv(Aoc{oci,1})'/(Aoc{oci,1})*[(Xk_temp(nn)-Xoc(oci));(Yk_temp(nn)-Yoc(oci))];
            lc = Roc(oci)*[Xk_temp(nn)-Xoc(oci);Yk_temp(nn)-Yoc(oci)]/norm( inv(Aoc{oci,1})*[Xk_temp(nn)-Xoc(oci);Yk_temp(nn)-Yoc(oci)]);
            lmove(nn,oci) = ( Omega0(nn,oci) + Dgc(:,nn,oci)'*(lc + [Xoc(oci);Yoc(oci)]-[Xk_temp(nn);Yk_temp(nn)]));
            end
        end
    end
    x  = sdpvar(1,N);
    y  = sdpvar(1,N);
    z  = sdpvar(1,N);
    vx = sdpvar(1,N-1);
    vy = sdpvar(1,N-1);
    vz = sdpvar(1,N-1);
    ux = sdpvar(1,N-2);
    uy = sdpvar(1,N-2);
    uz = sdpvar(1,N-2);
    T1  = sdpvar(1,1);
    T2  = sdpvar(1,1);
    F  = [  x(1,1) == x0];
    F  = F + [x(1,N) == xf];
    F  = F +  [y(1,1) == y0];
    F  = F +  [y(1,N) == yf];
    F  = F +  [z(1,1) == z0];
    F  = F +  [z(1,N) == zf];
    F  = F +  [vx(1,1) == Vx_0*T1];
    F  = F +  [vy(1,1) == Vy_0*T1];
    F  = F +  [vz(1,1) == Vz_0*T1];
    F  = F +  [vx(1,N-1) == Vx_f*T1];
    F  = F +  [vy(1,N-1) == Vy_f*T1];
    F  = F +  [vz(1,N-1) == Vz_f*T1];   
    %     F  = F +  [z >= 0];
    F  = F +  [.0 <= T1];
    F  = F +  [.0 <= T2];
    
    for i = 1:N-1
        % First-order Euler method
        F  = F + [x(1,i+1) == x(1,i) + tau*vx(1,i)];
        F  = F + [y(1,i+1) == y(1,i) + tau*vy(1,i)];
        F  = F + [z(1,i+1) == z(1,i) + tau*vz(1,i)];       
        F  = F + [Vzmin*T1 <= vz(1,i) <= Vzmax*T1];             % linear constraints for As/Descent speed
%         F  = F + [ norm([vx(1,i);vy(1,i)]) <= Vlmax*T1];        % cone constraints for level speed        
        F  = F + [ norm([vx(1,i);vy(1,i);vz(1,i)]) <= Speedmax*T1];  %second-order cone for velocity
%         F  = F + [ norm([vx(1,i);vy(1,i);vz(1,i)])^2 <= Vmax(1,i)*T2]; %rotate  cone for velocity
    end
    
    for j = 1:N-2
        F  = F + [vx(1,j+1) == vx(1,j) + tau*ux(1,j)];
        F  = F + [vy(1,j+1) == vy(1,j) + tau*uy(1,j)];
        F  = F + [vz(1,j+1) == vz(1,j) + tau*(uz(1,j) - g*T2/anorm)];
    end
    
    
    for ui = 1:N-2
        F  = F + [ norm([ux(1,ui);uy(1,ui);uz(1,ui)]) <= amax*T2];
        F  = F + [norm([ux(1,ui);uy(1,ui)]) <= tan(Tiltmax)*uz(1,ui)]; %Attitude angle constraint.
        F  = F + [T2*uzmin <= uz(1,ui)]; %Approximate minimum acceleration constraint.        
    end
    
    %     F  = F +  [ T2 >= T1^2 ];
    F  = F + rcone(T1,T2,0.5);
    
    %     F  = F +  [T2 == T1k_temp^2 + 2*T1k_temp*(T1-T1k_temp)];
    %     F  = F +  [T2 <= T1k_temp^2 + 2*T1k_temp*(T1-T1k_temp)];
    if flagobs == 1
        for ocj = 1:length(Roc)
            for oi = 1:N
                switch method
                    case 'Linearization'
                        F  = [F; Omega0(oi,ocj)+Dgc(:,oi,ocj)'*[(x(oi)-Xk_temp(oi));(y(oi)-Yk_temp(oi))]>= 0];
                    case 'Projection'
                        F  = [F; Omega0(oi,ocj)-lmove(oi,ocj)+Dgc(:,oi,ocj)'*[(x(oi)-Xk_temp(oi));(y(oi)-Yk_temp(oi))]>= 0];
                        %                     F  = [F; gc0(oi,ocj)-pc(oi,ocj)+Dgc(:,oi,ocj)'*[(x(oi)-Xk_temp(oi));(y(oi)-Yk_temp(oi))]>= 0];
                end
            end
        end
    end
    
    ops = sdpsettings('solver','mosek','verbose',1);%mosek ipopt ecos
    %     obj = T2 - Lambda(s)*T1;
    obj = T2 - Lambda(s)*T1;   
    %     obj = T2 - 0.5*Lambda(s)*T1;
%     obj = T2 - 0.3*14.4749386260*T1;     
    %     if s >= 12
    %         obj = T2 - Lambda(6)*T1;
    %     else
    %         obj = T2 - 0.3*Lambda(s)*T1;
    %     end
    solvesdp(F, obj, ops)
    TT2  = double(T2);
    TT1  = double(T1);
    actierror = TT1^2-TT2
    X   = double(x);
    Y   = double(y);
    Z   = double(z);
    Vx   = double(vx);
    Vy   = double(vy);
    Vz   = double(vz);
    Ux   = double(ux);
    Uy   = double(uy);
    Uz   = double(uz);
    Fobj(s) = double(obj);

    Xk(s+1,1:N) = X; Yk(s+1,1:N) = Y; Zk(s+1,1:N) = Z;
    vxk(s+1,:) = Vx; vyk(s+1,:) = Vy; vzk(s+1,:) = Vz;
    uxk(s+1,:) = Ux; uyk(s+1,:) = Uy; uzk(s+1,:) = Uz;
    T2k(s+1,1) = TT2;T1k(s+1,1) = TT1;
    errorx(s) = max(abs(Xk(s+1,:)-Xk(s,:)));
    errory(s) = max(abs(Yk(s+1,:)-Yk(s,:)));
    errorz(s) = max(abs(Zk(s+1,:)-Zk(s,:)));
    errort(s) = abs(T1k(s+1,:)-T1k(s,:));
    if errorx(s) <= tol && errory(s) <= tol && errorz(s) <= tol || flagobs == 0
        q_a = 1;
        q_b = -Lambda(s);
        q_c = Lambda(s)*TT1-TT2;
        T1_inters = 0.5*(-q_b+sqrt(q_b^2-4*q_a*q_c))/q_a;
        Lambda(s+1,:)= 2*T1_inters;%2*sqrt(T1_inters)
    else
        Lambda(s+1,:)=Lambda(s,:);
    end
    fprintf('SCOP completed %d times optimization,the optimal flight time is %d seconds.\n ',s,sqrt(TT2));
    fprintf('the errorx, errory and errorz are %d m, %d m and %d m.\n',errorx(s)*Lnorm,errory(s)*Lnorm,errorz(s)*Lnorm);
    fprintf('the errort is %d s.\n',errort(s)*Tnorm);
    %%%
    figure(11)
    plot3(Xk(end,:),Yk(end,:),Zk(end,:),'--b','Linewidth',1);
    xlabel('x(m)');ylabel('y(m)');zlabel('z(m)');
    hold on;grid on;
    time =linspace(0,tau*N,N)*TT1*Tnorm;
    Vv = sqrt(vxk(end,:).^2 + vyk(end,:).^2 + vzk(end,:).^2)/TT1;
    Vlev = sqrt(vxk(end,:).^2 + vyk(end,:).^2 )/TT1;
    figure(12)
    subplot(2,1,1)
    plot(time(1:end-1),Vv*Vnorm,'Linewidth',2);
    xlabel('time(s)');ylabel('V(m/s)');hold on;
    subplot(2,1,2)
    plot(time(1:end-1),Vlev*Vnorm,'Linewidth',2);
    xlabel('time(s)');ylabel('Vlev(m/s)');hold on;
    figure(22)
    subplot(3,1,1)
    plot(time(1:end-1),Vx*Vnorm,'Linewidth',2);
    xlabel('time(s)');ylabel('V_x(m/s)');hold on;
    subplot(3,1,2)
    plot(time(1:end-1),Vy*Vnorm,'Linewidth',2);
    xlabel('time(s)');ylabel('V_y(m/s)');hold on;
    subplot(3,1,3)
    plot(time(1:end-1),Vz*Vnorm,'Linewidth',2);
    xlabel('time(s)');ylabel('V_z(m/s)');hold on;
    if errorx(s) <= tol && errory(s) <= tol && errorz(s) <= tol...
            && errort(s) <= tol && abs(actierror) <=tol && flagobs == 1
        break
    elseif abs(actierror) <=tol && flagobs == 0
        break
    end
end
%%
for i = 1:N-2
    ux_norm=(Ux/TT2*anorm)./norm([Ux/(TT2)*anorm,Uy/(TT2)*anorm,Uz/(TT2)*anorm],2);
    uy_norm=(Uy/TT2*anorm)./norm([Ux/(TT2)*anorm,Uy/(TT2)*anorm,Uz/(TT2)*anorm],2);
    uz_norm=(Uz/TT2*anorm)./norm([Ux/(TT2)*anorm,Uy/(TT2)*anorm,Uz/(TT2)*anorm],2);
end
figure(1)
siz = size(Xk);
for i = 1:siz(1)-1
    plot3(Xk(i,:)*Lnorm,Yk(i,:)*Lnorm,Zk(i,:)*Lnorm,'--b','Linewidth',1);
    hold on;grid on;
end
plot3(Xk(end,:)*Lnorm, Yk(end,:)*Lnorm,Zk(end,:)*Lnorm,'r','Linewidth',2);
hold on;
if flagobs == 1
    ObstaclePlot;
end
%%
hold on;
figure(2)
labt = round(linspace(0,TT1*Tnorm,10));
labx = interp1(time,Xk(end,:)*Lnorm,labt,'linear','extrap');
laby = interp1(time,Yk(end,:)*Lnorm,labt,'linear','extrap');
labz = interp1(time,Zk(end,:)*Lnorm,labt,'linear','extrap');
plot3(Xk(end,:)*Lnorm,Yk(end,:)*Lnorm,Zk(end,:)*Lnorm,'r','Linewidth',2);hold on;
quiver3(Xk(end,3:end)*Lnorm,Yk(end,3:end)*Lnorm,Zk(end,3:end)*Lnorm,ux_norm,uy_norm,uz_norm,0.9,'b','Linewidth',0.5);
text(Xk(end,end)*Lnorm,Yk(end,end)*Lnorm,Zk(end,end)*Lnorm,num2str(TT1*Tnorm));
for i = 1:length(labt)-1
    text(labx(i),laby(i),labz(i),num2str(labt(i)))
end
xlim([min(Xk(end,:)*Lnorm) max(Xk(end,:)*Lnorm)]*1.1);
ylim([min(Yk(end,:)*Lnorm) max(Yk(end,:)*Lnorm)]*1.1);
axis equal;
grid on;
%%
figure(3)
subplot(3,1,1)
plot(time(1:N-1),Vx(1:N-1)/TT1*Vnorm,'k','Linewidth',2);
xlabel('time(s)');ylabel('V_x');xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
subplot(3,1,2)
plot(time(1:N-1),Vy(1:N-1)/TT1*Vnorm,'k','Linewidth',2)
xlabel('time(s)');ylabel('V_y');xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
subplot(3,1,3)
plot(time(1:N-1),Vz(1:N-1)/TT1*Vnorm,'k','Linewidth',2);
xlabel('time(s)');ylabel('V_z');xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
title([method ' ' 'Method']);
figure(4)
Vp=sqrt(Vx.^2+Vy.^2+Vz.^2)*Vnorm/TT1;
plot(time(1:N-1),Speedmax*Vnorm*ones(1,N-1),'r','Linewidth',2.5);
hold on;
plot(time(1:N-1),Vp(1:N-1),'--','Linewidth',1.5);
xlim([0,tau*N*TT1*Tnorm*1.05]);ylim([0,max(Vp)*1.2]);
xlabel('time(s)');ylabel('V(m/s)');grid on;
legend('Vconst','Vcone');title([method ' ' 'Method']);
phi = atan(Vz./Vy);
psi = asin(Vx./sqrt(Vx.^2+Vy.^2+Vz.^2));
figure(5)
subplot(2,1,1);
plot(time((1:N-1)),phi(1:N-1)*180/pi,'k','Linewidth',2);
xlabel('time(s)');ylabel('\phi');xlim([0,tau*N*TT1*Tnorm*1.05]);
grid on;title([method ' ' 'Method']);
subplot(2,1,2);
plot(time((1:N-2)),(diff(phi)/(TT1*Tnorm*tau))*180/pi,'k','Linewidth',2);
xlabel('time(s)');ylabel({'$\dot{\phi}$'},'Interpreter','latex');
xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
figure(6)
subplot(2,1,1);
plot(time((1:N-1)),psi(1:N-1)*180/pi,'k','Linewidth',2);
xlabel('time(s)');ylabel('\psi');xlim([0,tau*N*TT1*Tnorm*1.05]);
grid on;title([method ' ' 'Method']);
subplot(2,1,2);
plot(time((1:N-2)),(diff(psi)/(TT1*Tnorm*tau))*180/pi,'k','Linewidth',2);
xlabel('time(s)');ylabel({'$\dot{\psi}$'},'Interpreter','latex');
xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
% %% %
% figure(7)
% subplot(2,1,1)
% plot(time(1:N)*obj,Gaps(1:N),'k','Linewidth',2)
% xlabel('time(s)');ylabel({'$gap_{s}$'},'Interpreter','latex');
% xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
% subplot(2,1,2)
% plot(time(1:N)*obj,Gapc(1:N),'k','Linewidth',2)
% xlabel('time(s)');ylabel({'$gap_{c}$'},'Interpreter','latex');
% xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
% %%
% figure(8)
% subplot(3,1,1);
% plot(time,Xk(end,:)*Lnorm,'k','Linewidth',1.5);
% xlabel('time(s)');ylabel({'$x$'},'Interpreter','latex');
% xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
% subplot(3,1,2);
% plot(time,Yk(end,:)*Lnorm,'k','Linewidth',1.5);
% xlabel('time(s)');ylabel({'$y$'},'Interpreter','latex');
% xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
% subplot(3,1,3);
% plot(time,Zk(end,:)*Lnorm,'k','Linewidth',1.5);
% xlabel('time(s)');ylabel({'$z$'},'Interpreter','latex');
% xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
%%

% figure(9)
% plot3(Vx/TT1*Vnorm,Vy/TT1*Vnorm,Vz/TT1*Vnorm,'k','Linewidth',2);
% xlabel({'$V_x$'},'Interpreter','latex');
% ylabel({'$V_y$'},'Interpreter','latex');
% zlabel({'$V_z$'},'Interpreter','latex');
% grid on;axis equal;
hold on
figure(10)
subplot(2,1,1)
an = (diff(psi)/(TT1*Tnorm*tau)).*Vp(1:N-2);
at = (diff(Vp)/(TT1*Tnorm*tau));
plot(time((1:N-2)),sqrt(an.^2+at.^2),'r','Linewidth',2);
hold on;
plot(time(1:end-2),sqrt((Ux/TT2*anorm).^2+(Uy/TT2*anorm-1).^2+(Uz/TT2*anorm-1).^2),'--b','Linewidth',2);
grid on;
xlabel('time(s)');ylabel({'$a(t)$'},'Interpreter','latex');
legend('Speed change rate','True acceleration','Location','SouthEast');
subplot(2,1,2)
plot([time(1) ,time(end-2)],[Accel_max Accel_max],'-r','Linewidth',2);
hold on
plot(time(1:end-2),sqrt((Ux/TT2).^2+(Uy/(TT2)).^2+(Uz/(TT2)).^2)*anorm,'--b','Linewidth',2);
ylim([0,Accel_max*1.2]);grid on;
legend('Max Limit','Control','Location','SouthEast');
%%
hold on
figure(111)
t1x = (0:0.01:sqrt(max(T2k))+0.01)*Tnorm;
t2y = t1x.^2;
plot(t1x,t2y,'-r','Linewidth',2);
xlabel({'$t_1 (s)$'},'Interpreter','latex');
ylabel({'$t_2 (s^2)$'},'Interpreter','latex');hold on;grid on;
plot(T1k(2:end)*Tnorm,T2k(2:end)*Tnorm^2,'.--b','Linewidth',2,'MarkerSize',20);
figure(112);
subplot(2,2,1)
semilogy(errorx*Lnorm,'.-k','Linewidth',1.5,'MarkerSize',16);
grid on;xlim([0 s]);set(gca,'XTick',0:1:s);
xlabel('Iterations');ylabel({'$\delta_x (m)$'},'Interpreter','latex');
subplot(2,2,2);
semilogy(errory*Lnorm,'.-k','Linewidth',1.5,'MarkerSize',16);
grid on;xlim([0 s]);set(gca,'XTick',0:1:s);
xlabel('Iterations');ylabel({'$\delta_y (m)$'},'Interpreter','latex');
subplot(2,2,3);
semilogy(errorz*Lnorm,'.-k','Linewidth',1.5,'MarkerSize',16);
grid on;xlim([0 s]);set(gca,'XTick',0:1:s);
xlabel('Iterations');ylabel({'$\delta_z (m)$'},'Interpreter','latex');
subplot(2,2,4);
semilogy(errort*Tnorm,'.-k','Linewidth',1.5,'MarkerSize',16);
grid on;xlim([0 s]);set(gca,'XTick',0:1:s);
xlabel('Iterations');ylabel({'$\delta_t (t)$'},'Interpreter','latex');
%%
gamma = zeros(N-2,1);
for i = 1 :N-2
    U= [Ux(i),Uy(i),Uz(i)];
    V = [Vx(i),Vy(i),Vz(i)];
gamma(i) = acos(dot(U,V)/(norm(U)*norm(V)));
end
figure(324)
plot(time(1:end-2),gamma*180/pi,'--b','Linewidth',2);