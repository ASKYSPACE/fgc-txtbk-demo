figure(3131)
for i = 2:step
    Txk   = mk(i,2:end).*uxk(i,:);
    Tyk   = mk(i,2:end).*uyk(i,:);
    Tzk   = mk(i,2:end).*uzk(i,:);
    Tk =sqrt(Txk.^2 + Tyk.^2 + Tzk.^2);
plot((0:N-2)*tau,Tk,'--k', 'LineWidth', lineW*.5);
if i == step
plot((0:N-2)*tau,Tk,'k', 'LineWidth', lineW);
end
xlabel('time (s)');ylabel('throttle level (N)');grid on;hold on;
end