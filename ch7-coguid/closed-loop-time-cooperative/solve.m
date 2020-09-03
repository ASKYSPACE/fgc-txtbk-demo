function Y = solve(X)
global  r v step Tgo n N T0 h R0
t = h*step-2*h; % 当前时刻
tgo_d = 100-t;  % 期望的剩余攻击时间
nm = 5; %导弹最大可用过载
r = zeros(1,n);  % 弹目距离
sita = zeros(1,n); % 弹道倾角
q = zeros(1,n); % 视线角
sgm = zeros(1,n); % 前置角
for i = 1:n
r(i) = X(3*i);       %弹目距离
sita(i) = X(3*i-1); %导弹弹道倾角
q(i) = X(3*i-2);    %视线角
sgm(i) = sita(i) - q(i); %前置角
end
T_tt = 0; % 总剩余攻击时间
for i = 1:n
    Tgo(step,i)=( r(i)/v(i) )*( 1 + sgm(i)^2/( 4*N - 2 ) );
    T_tt = Tgo(step,i)+T_tt;
end
Tgo_ave = T_tt/n; %平均剩余攻击时间
% Tgo_ave = tgo_d; %期望攻击时间
FC = 0; % 剩余攻击时间方差
for i = 1:n
    FC = FC + (Tgo_ave-Tgo(step,i))^2;
end
BZC = sqrt(FC);  % 剩余攻击时间标准差
q1 = zeros(1,n);   %初始化视线角速度
sita1 = zeros(1,n);% 初始化弹道倾角速度
r1 = zeros(1,n);   % 初始化弹目逼近速度
for i = 1:n
    q1(i) = ( -v(i)*sin( sgm(i) ) - 0 )/r(i);  %视线角速度


    Ke = 40/( R0*T0 );  % Ke在前置角为0时会产生奇异
    omiga = Ke*r(i)*n*( Tgo_ave - Tgo(step,i) );
    


sita1(i)=N*( 1-omiga )*v(i)*q1(i);
sita1(i) = sita1(i)/v(i);
if (v(i)*sita1(i)/9.8)>nm
    sita1(i) = nm*9.8/v(i);
else if (v(i)*sita1(i)/9.8)< -nm
    sita1(i) = -nm*9.8/v(i);
    end
end
r1(i) = 0 - v(i)*cos( sgm(i) ); %求弹目距离导数
end
for i = 1:n
    if X(3*i)<1
        Y(3*i-2:3*i) = [0 0 0];
    else
        Y(3*i-2:3*i) = [ q1(i) sita1(i) r1(i) ];
    end
end







