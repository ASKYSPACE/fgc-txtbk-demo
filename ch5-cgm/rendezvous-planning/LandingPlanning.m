clc
clear
close all
cant_angle = 27;            % unit: deg
n = 6;                      % the number of identical thrusters
thrust_max_each = 3.1*10^3; % maximum thrust is 3.1kN for each thruster
Isp = 225;
g0 = 9.807;                 % Earth's gravitational constant
gm = 3.72;                  % Mars' gravitational constant 
m_wet = 1905;
m_dry = 1505;
alpha = 1/(Isp*g0*cos(cant_angle*pi/180)); % a constant for a given engine: m_dot = -c*||T||
thrust_max = n*thrust_max_each*cos(cant_angle*pi/180);
T_max = 0.8*thrust_max;                    % the maximum allowed thrust
T_min = 0.3*thrust_max;                    % the minimum allowed thrust
%% >>>>>>>>>> initial conditions <<<<<<<<< %%
x0  = 1.5e3; 
y0  = 0; 
z0  = 2.0e3;
Vx0 = -75; 
Vy0 = 0;  
Vz0 = 100;
%% >>>>>>>>>> terminal conditions <<<<<<<<< %%
tf = 75; % total flight time
xf  = 0;
yf  = 0; 
zf  = 0;
Vxf = 0; 
Vyf = 0;  
Vzf = 0;

N = 100;    % the number of discretized points
tau = tf/N; % time interval

m0 = m_wet;
lm0 = log(m0);
lmt0 = zeros(N,1);
for i=1:N
    lmt0(i) = log(m0-alpha*T_max*i*tau);
end
miu1 = T_min*exp(-lmt0);
miu2 = T_max*exp(-lmt0);

%% >>>>>>>>>>>>>±ﬂΩÁ”Î‘º ¯<<<<<<<<<<< %%

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
F  = [  x(1,1) == x0];
F  = F + [x(1,N) == xf];
F  = F +  [y(1,1) == y0];
F  = F +  [y(1,N) == yf];
F  = F +  [z(1,1) == z0];
F  = F +  [z(1,N) == zf];
% F  = F +  [x>=0];
F  = F +  [lm(1,1) == lm0];
F  = F +  [vx(1,1) == Vx0];
F  = F +  [vx(1,N) == Vxf];
F  = F +  [vy(1,1) == Vy0];
F  = F +  [vy(1,N) == Vyf];
F  = F +  [vz(1,1) == Vz0];
F  = F +  [vz(1,N) == Vzf];

for i = 1:N-1
    % First-order Euler method
    F  = F + [x(1,i+1)==x(1,i)+tau*vx(1,i)];
    F  = F + [y(1,i+1)==y(1,i)+tau*vy(1,i)];
    F  = F + [z(1,i+1)==z(1,i)+tau*vz(1,i)];
    F  = F + [lm(1,i+1)==lm(1,i)-tau*alpha*sig(1,i)];
    F  = F + [vx(1,i+1)==vx(1,i)+tau*(-0 + ux(1,i))];
    F  = F + [vy(1,i+1)==vy(1,i)+tau*uy(1,i)];
    F  = F + [vz(1,i+1)==vz(1,i)+tau*uz(1,i)]; 
end
for j = 1:N-1
    F  = F + [norm([ux(1,j);uy(1,j);uz(1,j)],2)<=sig(1,j)];    
    F  = F + [miu1(j)*(1-lm(1,j)+lmt0(j)+0.5*(lm(1,j)-lmt0(j))^2)<=sig(1,j)];
    F  = F + [sig(1,j)<=miu2(j)*(1-lm(1,j)+lmt0(j))];
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
Sigma= double(sig);
m    = exp(Lm);
Tx   = m(1:end-1).*Ux;
Ty   = m(1:end-1).*Uy;
Tz   = m(1:end-1).*Uz;
%% >>>>>>> plot the results <<<<<<<<< %%   
lineW = 2; 
figure(1)
subplot(2,1,1);
plot((0:N-1)*tau,X,'-r', 'LineWidth', lineW); grid on; hold on;
plot((0:N-1)*tau,Y,'-b', 'LineWidth', lineW); hold on;
plot((0:N-1)*tau,Z,'-k', 'LineWidth', lineW); hold on;
xlabel('time (s)'); ylabel('position (m)');
subplot(2,1,2);
plot((0:N-1)*tau,Vx,'-r', 'LineWidth', lineW); grid on; hold on;
plot((0:N-1)*tau,Vy,'-b', 'LineWidth', lineW); hold on;
plot((0:N-1)*tau,Vz,'-k', 'LineWidth', lineW); hold on;
xlabel('time (s)'); ylabel('velocity (m/s)');
figure(2)
plot3(X,Y,Z,'b', 'LineWidth', lineW);
xlabel('x(m)');ylabel('y(m)');zlabel('z(m)');
hold on;grid on;
Vv = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
figure(12)
subplot(3,1,1)
plot((0:N-1)*tau,Vx,'Linewidth',2);
xlabel('time(s)');ylabel('Vx(m/s)');grid on;hold on;
subplot(3,1,2)
plot((0:N-1)*tau,Vy,'Linewidth',2);
xlabel('time(s)');ylabel('Vy(m/s)');grid on;hold on;
subplot(3,1,3)
plot((0:N-1)*tau,Vz,'Linewidth',2);
xlabel('time(s)');ylabel('Vz(m/s)');grid on;hold on;
figure(13)
plot((0:N-1)*tau,Vv,'Linewidth',2);
xlabel('time(s)');ylabel('V(m/s)');grid on;hold on;
figure(3);
subplot(2,1,1);
plot((0:N-2)*tau,sqrt(Tx.^2+Ty.^2+Tz.^2)/1000,'-r', 'LineWidth', lineW);
xlabel('time (s)'); ylabel('throttle level'); grid on; 
subplot(2,1,2);
plot((0:N-1)*tau,m, 'LineWidth', lineW); 
xlabel('time (s)'); ylabel('mass (kg)'); grid on; 