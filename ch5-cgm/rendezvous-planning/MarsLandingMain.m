% dynamic equation of mass is considered in this program and T/m is viewed
% as the control
% Here are the data for your implementation (position in m, velocity in m/s)
% Initial condition
%% ************************code by Xinfu Liu (2015)********************* %%
%% ****** Please refer the paper by Behcet Acikmese: ******* %%
%% ****** "Convex Programming Approach to Powered Descent Guidance for Mars Landing", JGCD (2007) ******** %%

clear all;
close all;
format long;
cant_angle = 27; % unit: deg
n = 6; % the number of identical thrusters
thrust_max_each = 3.1*10^3; % maximum thrust is 3.1kN for each thruster
Isp = 225;
g0 = 9.807;                 % Earth's gravitational constant
gm = 3.72;                  % Mars' gravitational constant
m_wet = 1905;
m_dry = 1505;
alpha = 1/(Isp*g0*cos(cant_angle*pi/180)); % a constant for a given engine: m_dot = -c*||T||
thrust_max = n*thrust_max_each*cos(cant_angle*pi/180);
T_max = 0.8*thrust_max;  % the maximum allowed thrust
T_min = 0.3*thrust_max;  % the minimum allowed thrust
% T_min = 0;

nx = 7; % the number of elements in the state vector
nu = 4; % the number of elements in the control vector

%% ===== initial conditions ===== %%
r01 = 1.5*10^3;
r02 = 0;
r03 = 2*10^3;
v01 = -75;
v02 = 0;
v03 = 100;

%% ===== terminal conditions ===== %%
tf = 75; % total flight time
rf1 = 0;
rf2 = 0;
rf3 = 0;
vf1 = 0;
vf2 = 0;
vf3 = 0;

N = 100;  % the number of discretized points
Tau = tf/N; % time interval

AA = [zeros(3) eye(3) zeros(3,1);zeros(3) zeros(3) zeros(3,1);zeros(1,3) zeros(1,3) zeros(1)];
BB = [zeros(3) zeros(3,1);eye(3) zeros(3,1);zeros(1,3) -alpha];
bb = [0 0 0  -gm 0 0  0]';
H = eye(nx) + Tau*AA;
G = Tau*BB;

% mass = zeros(1,N);
% for i=1:N
%    mass(i) = m0-(m0-m_min)/tf*(i-1)*Tau
% end

%% ===== set an initial mass profile ====== %%
m0 = m_wet;
z0 = log(m0);
zz = zeros(N+1,1);
for i=1:N+1
    zz(i) = log(m0-alpha*T_max*(i-1)*Tau);
%    zz(i) = log(m_ini(i));
%    zz(i) = log(mass(i));
end
miu1 = T_min*exp(-zz);
miu2 = T_max*exp(-zz);

%% =====
ele = (N+1)*nx+N*nu; %x0,x1...xN, u0,u1...uN-1

%% ====== compute the necessary coefficients from discretization of the dynamics ==== %%
A = zeros(ele); 
A(1:nx,1:nx) = eye(nx);
b = zeros(ele,1);
b(1:nx) = zeros(nx,1);
for i=1:N
    b(i*nx+1:(i+1)*nx) = -Tau*bb; 
end
for i=1:N
   A(i*nx+1:(i+1)*nx,(i-1)*nx+1:i*nx) = H;  %coefficient of x(k-1)
   % row: nx + nx*(i-1)+1; column: nx*(i-1)+1
   A(i*nx+1:(i+1)*nx,(N+1)*nx+(i-1)*nu+1:(N+1)*nx+i*nu) = G; %coefficient of u(k-1)
   % row: nx + nx*(i-1)+1; column: (N+1)*nx+(i-1)*nu+1
   A((N+1)*nx+(i-1)*nu+1:(N+1)*nx+i*nu,(N+1)*nx+(i-1)*nu+1:(N+1)*nx+i*nu) = eye(nu);
   % row: (N+1)*nx+(i-1)*nu+1
%        A((N+1)*nx+N*nu+i,(N+1)*nx+N*nu+i) = 1;
end

%% ===== set the solver as mosek ===== %%
ops = sdpsettings('solver','ecos','verbose',1);

%% ===== define the optimization variables  ===== %%
y = sdpvar(ele,1); 

%% ===== define the objective function c*y ===== %%
c = zeros(1,ele);
c((N+1)*nx+nu:nu:ele) = 1;

%% ====== include all the constraints in F ======= %%
B = A-eye(ele);
F = [ B*y == b];  % dynamics
F = F + [y(1:nx) == [r01 r02 r03 v01 v02 v03 z0]', y(N*nx+1:N*nx+nx-1) == [rf1 rf2 rf3 vf1 vf2 vf3]']; % initial and terminal conditions
for i=1:N
  F = F + [norm( y((N+1)*nx+(i-1)*nu+1:(N+1)*nx+i*nu-1),2 ) <= y((N+1)*nx+i*nu)]; % norm of control           
  F = F + [y((N+1)*nx+i*nu) >= miu1(i)*( 1-(y(i*nx)-zz(i))+1/2*(y(i*nx)-zz(i))^2 )];
  F = F + [y((N+1)*nx+i*nu) <= miu2(i)*( 1-(y(i*nx)-zz(i)) )];
%   F = F + [ y(i*nx) >= log(m0-alpha*T_max*(i-1)*Tau)];
%   F = F + [ y(i*nx) <= log(m0-alpha*T_min*(i-1)*Tau)];
  F = F + [y((i-1)*nx+1) >= 0];  % this constraint is to make sure that the vehicle does not touch the ground
end         

%% ====== with constraint F and the objective function c*y, solve the problem ====== %%
solvesdp(F, c*y, ops);  % solve the optimization problem to get the solution y
y = double(y);
    
%% ====== plot the results ====== %%    
lineW = 2;    
figure(1)
subplot(2,1,1);
plot((0:N)*Tau,y(1:nx:(N+1)*nx),'-r', 'LineWidth', lineW); grid on; hold on;
plot((0:N)*Tau,y(2:nx:(N+1)*nx),'-b', 'LineWidth', lineW); hold on;
plot((0:N)*Tau,y(3:nx:(N+1)*nx),'-k', 'LineWidth', lineW); hold on;
xlabel('time (s)'); ylabel('position (m)');
subplot(2,1,2);
plot((0:N)*Tau,y(4:nx:(N+1)*nx),'-r', 'LineWidth', lineW); grid on; hold on;
plot((0:N)*Tau,y(5:nx:(N+1)*nx),'-b', 'LineWidth', lineW); hold on;
plot((0:N)*Tau,y(6:nx:(N+1)*nx),'-k', 'LineWidth', lineW); hold on;
xlabel('time (s)'); ylabel('velocity (m/s)');

figure(2);
subplot(2,1,1);
plot((0:N-1)*Tau,y((N+1)*nx+1:nu:ele),'-r', 'LineWidth', lineW); grid on; hold on;
plot((0:N-1)*Tau,y((N+1)*nx+2:nu:ele),'-b', 'LineWidth', lineW); hold on;
plot((0:N-1)*Tau,y((N+1)*nx+3:nu:ele),'-k', 'LineWidth', lineW); hold on;
xlabel('time (s)'); ylabel('acceleration (m/s^2)');
subplot(2,1,2);
mass = exp(y(nx:nx:N*nx));  % mass 
plot((0:N-1)*Tau,y((N+1)*nx+1:nu:ele).*mass/1000,'-r', 'LineWidth', lineW); grid on; hold on;
plot((0:N-1)*Tau,y((N+1)*nx+2:nu:ele).*mass/1000,'-b', 'LineWidth', lineW); hold on;
plot((0:N-1)*Tau,y((N+1)*nx+3:nu:ele).*mass/1000,'-k', 'LineWidth', lineW); hold on;
xlabel('time (s)'); ylabel('net force (kN)');

figure(3);
subplot(2,1,1);
plot((0:N-1)*Tau,y((N+1)*nx+nu:nu:(N+1)*nx+N*nu).*mass/thrust_max,'-r', 'LineWidth', lineW);
xlabel('time (s)'); ylabel('throttle level'); grid on; 
subplot(2,1,2);
plot((0:N)*Tau,exp(y(nx:nx:(N+1)*nx)), 'LineWidth', lineW); 
xlabel('time (s)'); ylabel('mass (kg)'); grid on; 
       
