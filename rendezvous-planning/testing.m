figure(6)
hold on;
axis([y0,yf,-40,30]);
drawsate = fill([-1.4,-1,-1,-1.4],[0.2,0.2,0.45,0.45],'g');
satellite = fill([-1.4,-1,-1,-1.4],[0.2,0.2,0.45,0.45],'g');
drawtext=text(0,0,'');
gaptext = {'12','21'};
timetext=text(0,0,'');
% fragment = [15 ,57 ,162 ,303]/tau;
% fragment = [15 ,45 ,84 ,303]/tau;

% fragment = [14 ,52 ,122 ,202]/tau;
fragment =[12 ,44 ,84 ,202]/tau; % off
frn = 1;
for i = 1:N
    for j = 1:length(sate_R)
        drawsate(j)=fill(sate_p(1,2,j)+sate_R(j)*cos(theta),sate_p(1,1,j)+sate_R(j)*sin(theta),...
            [1 0.618 1]);%[0.618 0.618 0.618]
        gaptext = num2str(round(distance(i,j),2));
        gaptext = strcat('Distance=',gaptext,'m');
%         drawtext(j)=text(sate_p(i,2,j)+6,sate_p(i,1,j)+6,gaptext);
%         set(drawtext(j),'Position',[sate_p(i,2,j)+sate_R(j),sate_p(i,1,j)-1.2*sate_R(j)]);
        
        set(drawsate(j),'EdgeColor', 'r');
        set(drawsate(j),'XData',sate_p(i,2,j)+sate_R(j)*cos(theta));
        set(drawsate(j),'YData',sate_p(i,1,j)+sate_R(j)*sin(theta));
        satellite(j)=fill(sate_p(1,2,j)+1*cos(theta),sate_p(1,1,j)+1*sin(theta),'b');
        set(satellite(j),'XData',sate_p(i,2,j)+1*cos(theta));
        set(satellite(j),'YData',sate_p(i,1,j)+1*sin(theta));
        
        draw1=plot(sate_p(1:i,2,j),sate_p(1:i,1,j),'-b', 'LineWidth', 0.75*lineW);hold on;
        delete(timetext);
        timetext = title(strcat('Timer:',num2str(round((i-1)*tau,2)),'s'),...
            'HorizontalAlignment','center');%'FontSize',12,
    end
    draw1=plot(Y(1:i),X(1:i),'k', 'LineWidth', lineW);
    xlabel('y(m)');ylabel('x(m)');grid on;
    set(gca,'XDir','reverse');
    pause(0.05);
    drawnow;
    F=getframe(gcf);
    if flagnote == 1
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,strcat(strtitle,'.gif'),'gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,strcat(strtitle,'.gif'),'WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;
    end
    if i==fragment(frn)
        saveas(gcf, strcat('T200obsoff_T',num2str((i-1)*tau)), 'fig');
        if  frn<length(fragment)
            frn = frn + 1;
        end
    end
    if i ~= N
        delete(satellite);
        delete(drawsate);
        delete(drawtext);
    end
    hold on;
end