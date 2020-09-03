%% 画剩余时间
t = (1:step)*h;
figure
for i = 1:n
plot(t,Tgo(:,i),'LineWidth',2);
hold on
end
title('导弹估计剩余攻击时间');
xlabel('时间/s');
ylabel('估计剩余时间/s');
grid on