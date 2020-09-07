function R=rk4_solve(f,a,b,Ya,N,M)
h=(b-a)/N;
T=zeros(1,N+1);
Y=zeros(M,N+1);
T=a:h:b;
Y(:,1)=Ya;
K1=zeros(M,1);
K2=zeros(M,1);
K3=zeros(M,1);
K4=zeros(M,1);


for j=1:N
    K1=feval(f,T(j),Y(:,j));
    K2=feval(f,T(j)+h/2,Y(:,j)+h*K1/2);
    K3=feval(f,T(j)+h/2,Y(:,j)+h*K2/2);
    K4=feval(f,T(j)+h,Y(:,j)+h*K3);
    Y(:,j+1)=Y(:,j)+h*(K1+2*K2+2*K3+K4)/6;
    
%     if Y(4,j+1)<0
%         break
%     end
end
% 
% for i=(N+1):-1:(j+1)
%         T(i)=[];
%         Y(:,i)=[];
%         
% end
        
R=[T' Y'];

%4阶龙格库塔原理与实现