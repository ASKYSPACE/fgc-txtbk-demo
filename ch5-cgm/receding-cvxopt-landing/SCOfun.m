%末段，使用P和beta，以高度y为自变量，成功，2018/7/4
function U_result=SCOfun(IV0,UC)
M=5; %最大迭代次数
N=70; %分段数

delta=10;


V0 = IV0(1);
ct0 = IV0(2);
x0 = IV0(3);
z0 = IV0(4);
m0 = IV0(5);


Vf=1;
ctf=-90*pi/180;
xf=529370+450;
zf=1221;
mf=59500;

V=linspace(V0,Vf,N);
ct=linspace(ct0,ctf,N);
x=linspace(x0,xf,N);
%z=linspace(z0,zf,N);
m=linspace(m0,mf,N);
P=linspace(1250000,1250000,N);
beta=linspace(0,0,N);

I=eye(1,N);
ERROR=[];

Isp=2958;
rou0=1.225;
S=pi*3.35*3.35/4;
b=0.0001208;
g=9.8;
nx=4;

H=(zf-z0)/(N-1);
z=z0:H:zf;
X0=[V0;ct0;x0;m0];
time=zeros(1,M);

% figure(1)
% plot(x,z);
% hold on;
%%
%tic

for j=1:M
    X = sdpvar(nx,N);
    U = sdpvar(2,N);  %P beta

    tic
    F = X(:,1)==X0;

    for i=1:N-1
        rou=rou0*exp(-b*z(i));
         ma=V(i)/340;

           Ma=[0.6 0.8 0.9 0.95 1.05 1.1 1.2 1.5 2 2.5 3 3.5 4 4.5 5];
           H1=0:20000:80000;
          zhi=[2.926	3.237	3.445	3.571	4.111	4.614	4.581	4.815	4.800	4.723	4.673	4.623	4.571	4.520	4.469
                3.080	3.407	3.626	3.759	4.327	4.857	4.822	5.069	5.053	4.972	4.919	4.866	4.812	4.758	4.704
                3.234	3.577	3.807	3.947	4.544	5.100	5.063	5.322	5.305	5.220	5.165	5.109	5.052	4.996	4.939
                3.696	4.088	4.351	4.511	5.193	5.828	5.786	6.083	6.063	5.966	5.903	5.839	5.774	5.709	5.645
                4.466	4.940	5.258	5.451	6.275	7.042	6.992	7.350	7.326	7.209	7.132	7.055	6.977	6.899	6.821];
          if ma<0.6
               ma=0.6;
          elseif ma>5
               ma=5;
         end
         Ca=interp2(Ma,H1,zhi,ma,z(i));
         
        D=0.5*rou*V(i)^2*S*Ca;
        
        f1=-(P(i)*cos(beta(i))+D)/(m(i)*V(i)*sin(ct(i)))-g/V(i);
        f2=-P(i)*sin(beta(i))/(m(i)*V(i)^2*sin(ct(i)))-g*cos(ct(i))/(V(i)^2*sin(ct(i)));
        f3=cos(ct(i))/sin(ct(i));
        f4=-P(i)/(Isp*V(i)*sin(ct(i)));        
        
        df1_dV=g/V(i)^2 + ((Ca*S*rou*V(i)^2)/2 + P(i)*cos(beta(i)))/(V(i)^2*m(i)*sin(ct(i))) - (Ca*S*rou)/(m(i)*sin(ct(i)));
        df1_dct=(cos(ct(i))*((Ca*S*rou*V(i)^2)/2 + P(i)*cos(beta(i))))/(V(i)*m(i)*sin(ct(i))^2);
        df1_dx=0;
        df1_dm=((Ca*S*rou*V(i)^2)/2 + P(i)*cos(beta(i)))/(V(i)*m(i)^2*sin(ct(i)));
        
        df2_dV=(2*g*cos(ct(i)))/(V(i)^3*sin(ct(i))) + (2*P(i)*sin(beta(i)))/(V(i)^3*m(i)*sin(ct(i)));
        df2_dct=g/V(i)^2 + (g*cos(ct(i))^2)/(V(i)^2*sin(ct(i))^2) + (P(i)*cos(ct(i))*sin(beta(i)))/(V(i)^2*m(i)*sin(ct(i))^2);
        df2_dx=0;
        df2_dm=(P(i)*sin(beta(i)))/(V(i)^2*m(i)^2*sin(ct(i)));
        
        df3_dV=0;
        df3_dct=- cos(ct(i))^2/sin(ct(i))^2 - 1;
        df3_dx=0;
        df3_dm=0;
        
        df4_dV=P(i)/(Isp*V(i)^2*sin(ct(i)));
        df4_dct=(P(i)*cos(ct(i)))/(Isp*V(i)*sin(ct(i))^2);
        df4_dx=0;
        df4_dm=0;
        
        df1_dP=-cos(beta(i))/(V(i)*m(i)*sin(ct(i)));
        df2_dP=-sin(beta(i))/(V(i)^2*m(i)*sin(ct(i)));
        df3_dP=0;
        df4_dP=-1/(Isp*V(i)*sin(ct(i)));
        
        df1_dbeta=(P(i)*sin(beta(i)))/(V(i)*m(i)*sin(ct(i)));
        df2_dbeta=-(P(i)*cos(beta(i)))/(V(i)^2*m(i)*sin(ct(i)));
        df3_dbeta=0;
        df4_dbeta=0;

        
        AA=[df1_dV, df1_dct, df1_dx, df1_dm
            df2_dV, df2_dct, df2_dx, df2_dm
            df3_dV, df3_dct, df3_dx, df3_dm
            df4_dV, df4_dct, df4_dx, df4_dm
            ];
        BB1=[df1_dP; df2_dP; df3_dP; df4_dP];
        BB2=[df1_dbeta; df2_dbeta; df3_dbeta; df4_dbeta];

        CC=[f1-df1_dV*V(i)-df1_dct*ct(i)-df1_dx*x(i)-df1_dm*m(i)-df1_dP*P(i)-df1_dbeta*beta(i);
            f2-df2_dV*V(i)-df2_dct*ct(i)-df2_dx*x(i)-df2_dm*m(i)-df2_dP*P(i)-df2_dbeta*beta(i);
            f3-df3_dV*V(i)-df3_dct*ct(i)-df3_dx*x(i)-df3_dm*m(i)-df3_dP*P(i)-df3_dbeta*beta(i);
            f4-df4_dV*V(i)-df4_dct*ct(i)-df4_dx*x(i)-df4_dm*m(i)-df4_dP*P(i)-df4_dbeta*beta(i);
            ];
        BB=[BB1,BB2];
        A=H*AA+eye(nx);
        B1=H*BB1;
        B2=H*BB2;

        C=H*CC;
        
        F = [F;X(:,i+1)== A*X(:,i) + B1*U(1,i)+B2*U(2,i)+ C];    %路径约束
        
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %总约束
                F = [F;0 <= X(1,:) <= V0];    %速度约束
        %        F = [F;-150*pi/180 <= X(2,:)<= -20*pi/180];   %弹道倾角约束
                F = [F;x0-100 <= X(3,:) <= xf+20];   %范围约束
%                F = [F;zf<= X(4,:) <=z0];       %高度约束
                F = [F;UC(1) <= U(1,1) <= UC(1)]; 
%                 F = [F;UC(2) <= U(2,1) <= UC(2)]; 

                %控制量约束
                for ii=1:N
                    F = [F;1250000 <= U(1,ii)<= 2500000];    %轴向推力约束
                    F = [F;-8*pi/180 <= U(2,ii) <= 8*pi/180];    %侧向推力约束

                    
                end
                for ii=1:N-1
                    F = [F;-0.01 <= X(2,ii+1)-X(2,ii) <= 0.01];
                    F = [F;-20000 <= U(1,ii+1)-U(1,ii) <= 20000];
                    F = [F;-0.01 <= U(2,ii+1)-U(2,ii) <= 0.01];
                end

                %末端约束
                F = [F; 0<=X(1,N)<=1];    %速度约束
                F = [F;58500<=X(4,N)];   %质量约束
                F = [F;-90*pi/180<=X(2,N)<=-90*pi/180];    %落角约束
                F = [F;xf<=X(3,N)<=xf];
               F = [F;U(2,N) ==0];        %推力方向末端
%                 F = [F;U(2,N)+X(2,N) == -90*pi/180];
    
    obj=-X(4,N);  %sum(U(3,:))

    time(j)=toc;
    warning('off','YALMIP:strict') ;
    ops = sdpsettings('solver','ecos','verbose',1);
    solvesdp(F, obj, ops)
    
    Xk(:,:,j) = double(X);
    Uk(:,:,j) = double(U);
    
%     figure(1);
%     plot(Xk(3,:,j),z);hold on; %弹道
    
    error=norm((abs(Xk(:,:,j)-[V;ct;x;m])),1)   %本次最优结果与上次最优结果对比 收敛判据 
    ERROR=[ERROR,error];
    if ((error < delta))  % && (Xk(3,end,j)<=xf+2) && (Xk(3,end,j)>=xf-2 )
        break;
    end
    %本次寻优结果赋给V ct x z m作为下次的离散矩阵的计算来源
    
        V=Xk(1,:,j);
        ct=Xk(2,:,j);
        x=Xk(3,:,j);
        m=Xk(4,:,j);
        P=Uk(1,:,j);
        beta=Uk(2,:,j);
         
end
%toc
X = double(X);
U = double(U);


%%
% t_elapsed=sum(sum(time))

t=zeros(1,N);
for iii=2:N
   %t(iii)=t(iii-1)+H/sin((X(2,iii)+X(2,iii-1))/2)/((X(1,iii)+X(1,iii-1))/2); 
   %t(iii)=t(iii-1)+H/sin(X(2,iii))/X(1,iii); 
   t(iii)=t(iii-1)+H/sin(X(2,iii-1))/X(1,iii-1); 
end
% figure(1);
% hold on;
% pp = plot(X(3,1:N),z,'r');

%ylim([0,6000]);
% xl = xlabel('Range /m');
% yl = ylabel('Hight /m');
% title('Optimized sequence trajectory');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16);
% set(pp,'LineWidth',2);
% grid on;
% hold on;



% figure(2);
% hold on;
% pp = plot(t,X(1,1:N),'r');
% 
% xl = xlabel('t /s');
% yl = ylabel('V /m/s');
% title('速度');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16);
% set(pp,'LineWidth',2);
% grid on;
% hold on;
% 
% 
% 
% figure(3);
% hold on;
% pp = plot(t,X(2,1:N)*180/pi,'r');
% 
% xl = xlabel('t /s');
% yl = ylabel('theta /°');
% title('弹道倾角');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16);
% set(pp,'LineWidth',2);
% grid on;
% hold on;
% 
% figure(4);
% hold on;
% pp = plot(t,X(4,1:N),'r');
% 
% xl = xlabel('t /s');
% yl = ylabel('m /kg');
% title('质量');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16);
% set(pp,'LineWidth',2);
% grid on;
% hold on;
% 
% figure(5);
% hold on;
% pp = plot(t,U(1,1:N)/1000,'r');
% 
% xl = xlabel('t /s');
% yl = ylabel('P /kN');
% title('推力大小');
% 
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16);
% set(pp,'LineWidth',2);
% grid on;
% hold on;
% 
% figure(6);
% hold on;
% pp = plot(t,U(2,1:N)*180/pi,'r');
% 
% xl = xlabel('t /s');
% yl = ylabel('alpha /°)');
% title('攻角');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16);
% set(pp,'LineWidth',2);
% grid on;
% hold on;
% 
% figure(9);
% hold on;
% pp = plot(X(3,1:N),z,'r');
% 
% %ylim([0,6000]);
% xl = xlabel('x /m');
% yl = ylabel('y /m');
% title('优化弹道');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16);
% set(pp,'LineWidth',2);
% grid on;
% hold on;
% 
% figure(10)
% plot(ERROR,'r');hold on;

%%

t_P=[U;t];
U_result = t_P;

% save('t_P.mat','t_P');
end
