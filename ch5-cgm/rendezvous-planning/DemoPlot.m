clc
% clear 
close all
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
    view([25,21]);
    grid on;xlabel('x(m)','FontName','Times New Roman','FontSize',11);
    ylabel('y(m)','FontName','Times New Roman','FontSize',11);
    zlabel('z(m)','FontName','Times New Roman','FontSize',11);
end
%%
hold on;
figure(2)
qurin = 1;
labt = round(linspace(0,TT1*Tnorm,12));
labx = interp1(time,Xk(end,:)*Lnorm,labt,'linear','extrap');
laby = interp1(time,Yk(end,:)*Lnorm,labt,'linear','extrap');
labz = interp1(time,Zk(end,:)*Lnorm,labt,'linear','extrap');
plot3(Xk(end,:)*Lnorm,Yk(end,:)*Lnorm,Zk(end,:)*Lnorm,'r','Linewidth',2);hold on;
quiver3(Xk(end,3:qurin:end)*Lnorm,Yk(end,3:qurin:end)*Lnorm,Zk(end,3:qurin:end)*Lnorm,...
    ux_norm(1:qurin:end),uy_norm(1:qurin:end),uz_norm(1:qurin:end),0.9,'b','Linewidth',0.5);
text(Xk(end,end)*Lnorm,Yk(end,end)*Lnorm,Zk(end,end)*Lnorm,[num2str(TT1*Tnorm),'s']...
    ,'FontName','Times New Roman','FontSize',11);
for i = 1:length(labt)-1
    text(labx(i),laby(i),labz(i),[num2str(labt(i)),'s']...
        ,'FontName','Times New Roman','FontSize',11);
end
xlim([min(Xk(end,:)*Lnorm) max(Xk(end,:)*Lnorm)]*1.1);
ylim([min(Yk(end,:)*Lnorm) max(Yk(end,:)*Lnorm)]*1.1);
axis equal;view([25,21]);
grid on;xlabel('x(m)','FontName','Times New Roman','FontSize',11);
ylabel('y(m)','FontName','Times New Roman','FontSize',11);
zlabel('z(m)','FontName','Times New Roman','FontSize',11);
%%
hold on
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
%%
figure(4)
Vp=sqrt(Vx.^2+Vy.^2+Vz.^2)*Vnorm/TT1;
plot(time(1:N-1),Speed_max*ones(1,N-1),'r','Linewidth',2.5);hold on;
plot(time(1:N-1),Vp(1:N-1),'--','Linewidth',1.5);
xlim([0,tau*N*TT1*Tnorm*1.05]);ylim([0,max(Vp)*1.2]);
xlabel('Time(s)','FontName','Times New Roman','FontSize',11);
ylabel('Velocity (m/s)','FontName','Times New Roman','FontSize',11);
grid on;
legend('The maximum velocity','Actual velocity');
% title([method ' ' 'Method']);
phi = atan(Vz./Vy);
psi = asin(Vx./sqrt(Vx.^2+Vy.^2+Vz.^2));
%%

figure(5)
hold on;
subplot(2,1,2);
plot(time((1:N-1)),phi(1:N-1)*180/pi,'k','Linewidth',2);
xlabel('Time(s)','FontName','Times New Roman','FontSize',11);
ylabel('Flight heading angle (deg)','FontName','Times New Roman','FontSize',11);
grid on;xlim([0,tau*N*TT1*Tnorm*1.05]);
hold on;
subplot(2,1,1);
plot(time((1:N-1)),psi(1:N-1)*180/pi,'k','Linewidth',2);
xlabel('Time(s)','FontName','Times New Roman','FontSize',11);
ylabel('Flight path angle (deg)','FontName','Times New Roman','FontSize',11);
grid on;xlim([0,tau*N*TT1*Tnorm*1.05]);
%%
figure(6)
subplot(2,1,1);
plot(time((1:N-2)),(diff(phi)/(TT1*Tnorm*tau))*180/pi,'k','Linewidth',2);
xlabel('time(s)');ylabel({'$\dot{\phi}$'},'Interpreter','latex');
xlim([0,tau*N*TT1*Tnorm*1.05]);grid on;
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
% subplot(2,1,1)
% an = (diff(psi)/(TT1*Tnorm*tau)).*Vp(1:N-2);
% at = (diff(Vp)/(TT1*Tnorm*tau));
% plot(time((1:N-2)),sqrt(an.^2+at.^2),'r','Linewidth',2);
% hold on;
% plot(time(1:end-2),sqrt((Ux/TT2*anorm).^2+(Uy/TT2*anorm-1).^2+(Uz/TT2*anorm-1).^2),'--b','Linewidth',2);
% grid on;
% xlabel('time(s)');ylabel({'$a(t)$'},'Interpreter','latex');
% legend('Speed change rate','True acceleration','Location','SouthEast');
% subplot(2,1,2)

plot(time(1:end-2),sqrt((Ux/TT2).^2+(Uy/(TT2)).^2+(Uz/(TT2)).^2)*anorm,'-b','Linewidth',2);
ylim([0,Accel_max*1.2]);grid on;hold on;
plot([time(1) ,time(end-2)],[Accel_max Accel_max],'--r','Linewidth',2);
legend('The maximum acceleration','Actual acceleration','Location','SouthEast');
xlabel('Time(s)','FontName','Times New Roman','FontSize',11);
ylabel('Accel. (m/s^2)','FontName','Times New Roman','FontSize',11);
%%
hold on
figure(111)
t1x = (0:0.01:sqrt(max(T2k))+0.01)*Tnorm;
t2y = t1x.^2;
plot(t1x,t2y,'-r','Linewidth',2);
xlabel({'$t_1 (s)$'},'Interpreter','latex');
ylabel({'$t_2 (s^2)$'},'Interpreter','latex');hold on;grid on;
plot(T1k(2:end)*Tnorm,T2k(2:end)*Tnorm^2,'.--b','Linewidth',2,'MarkerSize',20);
%%
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
figure(123)
PHI = acos(Uz./sqrt(Ux.^2+Uy.^2+Uz.^2));
plot(time(1:end-2),PHI*180/pi,'-b','Linewidth',2);grid on;
hold on;
plot([time(1) ,time(end-2)],[Tiltmax Tiltmax]*180/pi,'--r','Linewidth',2.5);
xlabel('Time(s)','FontName','Times New Roman','FontSize',11);
ylabel('Inclination angle (deg)','FontName','Times New Roman','FontSize',11);
ylim([0,Tiltmax*180/pi*1.05]);
% legend('The maximum inclination angle','Actual inclination angle');
legend('The maximum inclination angle','Actual inclination angle');
%%
gamma = zeros(N-2,1);
Ut   = zeros(N-2,1);
Un   = zeros(N-2,1);
Gt    = zeros(N-2,1);
Gn    = zeros(N-2,1);
for i = 1 :N-2
    U= [Ux(i),Uy(i),Uz(i)]/TT2;
    V = [Vx(i),Vy(i),Vz(i)/TT1];
gamma(i) = acos(dot(U,V)/(norm(U)*norm(V)));
Ut(i)    = norm(U)*cos(gamma(i))*anorm;
Un(i)    = norm(U)*sin(gamma(i))*anorm;
Gt(i)    = dot([0,0,-g],V/norm(V));
Gn(i)    =sqrt(g^2-Gt(i)^2);
end
figure(1236)
plot(time(1:end-2),gamma*180/pi,'-b','Linewidth',2);grid on;
xlabel('Time(s)','FontName','Times New Roman','FontSize',11);
ylabel('gamma angle (deg)','FontName','Times New Roman','FontSize',11);
figure(52)
subplot(2,1,1)
plot(time(1:end-2),Ut+Gt,'-b','Linewidth',2);grid on;
subplot(2,1,2)
plot(time(1:end-2),Un+Gn,'-b','Linewidth',2);grid on;