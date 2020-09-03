clear
clc
global h Pm Pt v n N step Tgo T0 R0 acc acc_b
n = 2; % 导弹数量
N = 3; % 比例导引系数
tf =160; % 最大仿真时间限制
acc = zeros(100,n); %
acc_b = zeros(100,n);

Pt = [14000,0];    %目标的位置
v = zeros(n,1);   %n个导弹的速度
Pm = zeros(n,2);  %导弹的位置(n个导弹)
Tgo = zeros(10,n); % 个导弹剩余攻击时间
q0 = zeros(1,n);  %初始化各视线角
Sgm = zeros(1,n);  %前置角
sita = zeros(1,n); %弹道倾角
r0 = zeros(1,n);  %初始距离
R_total = 0;
T_total = 0;
tgo = zeros(1,n);
for i = 1:n
    v(i) = 200+(i-1)*20;  % 导弹速度
    Pm(i,:) = [-(i-1)*1000,(i-1)*(00)]; % 导弹位置
    q0(i) = atan( ( Pt(2) - Pm(i,2) )/( Pt(1) - Pm(i,1) ) );%弹i视线角
    Sgm(i) = ( 30+(i-1)*15 )*pi/180;    %前置角
    sita(i) = q0(i)+Sgm(i);    %弹道倾角
    r0(i) = sqrt( ( Pt(2) - Pm(i,2) )^2 + ( Pt(1) - Pm(i,1) )^2 );%初始距离
    Tgo(1,i) = ( r0(i)/v(i) )*( 1 + ( Sgm(i)^2/( 2*( 2*N - 1 ) ) ) );  % 剩余攻击时间估计           
    R_total = R_total + r0(i);     %总距离
    T_total = T_total + Tgo(1,i);   % 总剩余攻击时间
    
end
R0 = R_total/n;   % 初始平均弹目距离
T0 = T_total/n;   % 初始平均剩余攻击时间
h = 0.005;      %仿真步长
t = 0:h:tf;     %时间
Step = 1+tf/h;  %最大步长
% 初始化状态量
X = zeros(10,3*n);   %状态量
for i = 1:n
    X(1,(3*i-2):(3*i)) = [ q0(i),sita(i),r0(i) ];
end
    % [ 弹i视线角，弹i弹道倾角，弹i弹目距离 ]
i = 1;
P = zeros( n,2,100 );  % 导弹的位置记录
flag = 1;  % 判断是否所有的导弹都命中目标
while (i <=Step)&&flag  %( norm( [ ( X(8)-X(6) ),( X(7)-X(5) ) ] )<10 )
    i = i+1;
    step = i;
    X(i,:) = RK_4( X(i-1,:) );
    for j = 1:n
        P(j,:,i) =  Pt - [ X(i,3*j)*cos( X(i,j*3-2) ),X(i,3*j)*sin( X(i,3*j-2) ) ];%记录导弹位置
    end
    flag = 1;
    for j = 1:n
        if X(i,3*j)<5
            flag = 1*flag;
        else
            flag = 0;
        end
        
    end
    if flag == 1
        flag = 0;
    else
        flag = 1;
    end
end
plotTgo
plotSgm
plotAcc
out
% plotdetaAcc


