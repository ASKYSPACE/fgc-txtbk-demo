% This matlab program is based on rolling convex optimization to
% solve the problem of rocket soft landing
% inputs: rocket velocity
%         rocket flight-angle 
%         rocket position
%         the mass of rocket
%         the time of guidance period
%         the number of guidance
%         the constraints of control
%         the constraints of terminal state
%         the constraints of process
clear all;
IV0 = [253,-68.45/57.3,529370,3031,66917]; %the initial state of rocket
T = 2;                     %the time of guidance period
M = 7;                     %the number of guidance
U_save = [];
t_now = 0;
UC = [1250000;1250000];
for i = 1:M
    U=SCOfun(IV0,UC);     
    P = U(1,:);
    alpha = U(2,:);
    t_end = U(3,end);
    U(3,:) = U(3,:)+t_now;
    save('U.mat','U');
    U_now = U;
    if i == M
        T = t_end;
    end
    
    while U_now(3,end)>t_now+T
       U_now(:,end) = [];
    end
    U_save = [U_save,U_now];  
    UC = U_now(1:2,end);
    
    R4=rk4_solve('jianmo',t_now,t_now+T,IV0,400,5);  
    t_now = t_now+T;
    IV0 = R4(end,2:end)
    t = R4(:,1);
    %
    figure(1);
    plot(R4(:,4),R4(:,5),'Color',[rand(),rand(),rand()],'LineWidth',1);
    hold on
    grid on
    
    figure(2)
    plot(t,R4(:,2),'Color',[rand(),rand(),rand()],'LineWidth',3);
    hold on
    grid on
%     
    figure(3)
    plot(t,R4(:,3)*180/pi,'Color',[rand(),rand(),rand()],'LineWidth',3);
    hold on
    grid on
    
    figure(4)
    plot(t,R4(:,6),'Color',[rand(),rand(),rand()],'LineWidth',3);
    hold on
    grid on
    
    figure(5)
    plot(U(3,:),U(1,:)/1000,'Color',[rand(),rand(),rand()],'LineWidth',3);
    hold on
    grid on
    
    figure(6)
    plot(U(3,:),U(2,:)*180/pi,'Color',[rand(),rand(),rand()],'LineWidth',3);
    hold on
    grid on
end
 figure(7)
 plot(U_save(3,:),U_save(1,:),'r','LineWidth',3);
 title('滚动推力曲线结果')
 hold on
 grid on
 figure(8)
 plot(U_save(3,:),U_save(2,:)*180/pi,'r','LineWidth',3);
 title('滚动攻角曲线结果')
 hold on
 grid on

