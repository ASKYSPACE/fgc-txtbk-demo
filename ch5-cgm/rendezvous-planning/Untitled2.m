clc
pic_num = 1;
I=imread('satellite2.jpg');
[xis,yis,zis]=size(I);
theta=0:pi/100:2*pi;
figure(5)
hold on;
axis([y0,yf,-40,30]);
drawsate = fill([-1.4,-1,-1,-1.4],[0.2,0.2,0.45,0.45],'g');
satellite = fill([-1.4,-1,-1,-1.4],[0.2,0.2,0.45,0.45],'g');
for i = 1:N
    draw1=plot(Y(1:i),X(1:i),'k', 'LineWidth', lineW);
    xlabel('y(m)');ylabel('x(m)');grid on;
    %     axis equal;
    set(gca,'XDir','reverse');
    %     axis image;
    for j = 1:length(sate_R)
        drawsate(j)=fill(sate_p(1,2,j)+sate_R(j)*cos(theta),sate_p(1,1,j)+sate_R(j)*sin(theta),...
            [1 0.618 1]);%[0.618 0.618 0.618]
        set(drawsate(j),'EdgeColor', 'r');
        set(drawsate(j),'XData',sate_p(i,2,j)+sate_R(j)*cos(theta)); 
        set(drawsate(j),'YData',sate_p(i,1,j)+sate_R(j)*sin(theta)); 
        satellite(j)=fill(sate_p(1,2,j)+1*cos(theta),sate_p(1,1,j)+1*sin(theta),'b');
        set(satellite(j),'XData',sate_p(i,2,j)+1*cos(theta)); 
        set(satellite(j),'YData',sate_p(i,1,j)+1*sin(theta)); 
        draw1=plot(sate_p(1:i,2,j),sate_p(1:i,1,j),'--b', 'LineWidth', lineW);hold on;  
    end
    pause(0.05);
    drawnow;
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,'test2.gif','gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,'test2.gif','gif','WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;
    if i ~= N
        delete(satellite);
        delete(drawsate);
    end
    hold on; 
end
%%

theta=0:pi/100:2*pi;
figure(6)
distance = zeros(N,length(sate_R));
for j = 1:length(sate_R)
    distance(:,j) = sqrt((Y'-sate_p(:,2,j)).^2 + (X'-sate_p(:,1,j)).^2)-sate_R(j);
    plot((0:N-1)*tau,distance(:,j),'k','LineWidth', lineW);
    hold on;grid on;
    xlabel('time (s)');ylabel('Distance (m)');
end
%   legend();
