% The function is used for trajectory integration. 
% For specified sliding mode variable structure guidance law, given initial
% state value and integration time, one can obtained trajectory. 

%  parameters and variables %
% c         coefficient of general switching function
% k         reaching law coefficient of general sliding-mode law
% epsi      switching term coefficient of general sliding-mode law
% k_1       coefficient of time-to-go switching function
% k_2       reaching law coefficient of time-to-go sliding-mode law
% kapa      switching term coefficient of time-to-go sliding-mode law
% lamda_d   desired line of sight angle
% V_0       initial velocity
% H_0       initial height
% x_0       initial horizontal distance
% theta_m0  flight path angles
% t0        integration time
%% 
% Guidance parameters
c = 0.2;
k = 0.2;
epsi = 0.002;
k_1 = 1;
k_2 = 1;
kapa = 0.01;
% Desired line of sight angle
lamda_d = -30*pi/180;
% Initial value
V_0 = 260;
H_0 = 3000;
x_0 = -6000;
theta_m0 = 0;
% Integration time
t0 = 28;
%% Trajectory simulation
% Calculate initial state value
r0 = sqrt(H_0.^2+x_0.^2);
lamda0 = atan(H_0./x_0);
theta_t0 = 0;
x_m0 = x_0;
H_m0 = H_0;
d_lamda0 = V_0*sin(lamda0-theta_m0)/r0;
s0 = d_lamda0+c*(lamda0-lamda_d);
y0 = [r0 lamda0 theta_t0 theta_m0 x_m0 H_m0 d_lamda0 s0];
t_go0 = r0/(V_0*cos(lamda0-theta_m0));
s_go0 = d_lamda0+k_1/t_go0*(lamda0-lamda_d);
y_go0 = [r0 lamda0 theta_t0 theta_m0 x_m0 H_m0 d_lamda0 s_go0];
% trajectory integration
h = 0.001;
t1=0:h:t0;
y_simple = ode4A(@DE1,t1,y0,V_0,c,k,epsi,lamda_d);
t_simple = t1(1:length(y_simple));
y_go = ode4A(@DE2,t1,y_go0,V_0,k_1,k_2,kapa,lamda_d);
t_go = t1(1:length(y_go));
%% Plot figures
% Note: you need to modify the parameters to recalculate.
r_go2 = y_go(:,1);r_simple2 = y_simple(:,1);
r_go1 = [r0;r_go2];r_simple1 = [r0;r_simple2];
r_go1(end) = [];r_simple1(end) = [];
dr_go = (r_go2-r_go1)/h;dr_simple = (r_simple2-r_simple1)/h;
tt_go = -r_go2./dr_go;
s_go = y_go(:,7)+k_1./tt_go.*(y_go(:,2)-lamda_d);
s_simple = y_simple(:,7)+c.*(y_simple(:,2)-lamda_d);
dr_go = (V_0.*cos(y_go(:,2)-y_go(:,4)));
A_M_go = ( y_go(:,1).*(-2*dr_go./y_go(:,1)+(k_1+k_2)./tt_go).*y_go(:,7)+k_1*y_go(:,1)./tt_go.^2*(k_2+1).*(y_go(:,2)-lamda_d)+...
    kapa*y_go(:,1)./tt_go.*sat(s_go) )./cos(y_go(:,2)-y_go(:,4));
A_M_simple = ( y_simple(:,1).*(-2*dr_simple./y_simple(:,1)+c+k).*y_simple(:,7)+c*y_simple(:,1)*k.*(y_simple(:,2)-lamda_d)+...
    epsi*y_simple(:,1).*sat(s_simple) )./cos(y_simple(:,2)-y_simple(:,4));

figure(1)
plot(y_go(:,5),y_go(:,6),'-')
title('弹道')
grid on;hold on
plot(y_simple(:,5),y_simple(:,6),'--')
legend('自适应','普通')
figure(2)
plot(t_go,y_go(:,2)*180/pi,'-')
grid on;hold on
plot(t_simple,y_simple(:,2)*180/pi,'--')
legend('自适应','普通')
title('视线角/。')
grid on

figure(3)
plot(t_go,A_M_go,'-')
title('加速度')
grid on;hold on
plot(t_simple,A_M_simple,'--')
legend('自适应','普通')

figure(4)
title('切换函数值')
grid on;hold on
plot(t_go,y_go(:,8),'-')
plot(t_simple,y_simple(:,8),'--')
legend('自适应','普通')