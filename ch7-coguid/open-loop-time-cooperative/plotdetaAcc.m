%% 附加加速度
figure
t = (1:step-1)*h;
for i = 1:n
plot(t,a(2:step,i),'LineWidth',2);
hold on
end
title('附加导弹加速度');
xlabel('时间/秒');
ylabel('加速度/(m/s^2)');
grid on