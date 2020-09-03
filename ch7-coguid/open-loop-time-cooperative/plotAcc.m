%% 画加速度

t = (1:step-1)*h;
for i = 1:n
figure
plot(t,acc(2:step,i),'LineWidth',2);
hold on
plot(t,acc_b(2:step,i),'LineWidth',2);
hold on
title('导弹加速度');
xlabel('时间/秒');
ylabel('加速度/(m/s^2)');
grid on
end

