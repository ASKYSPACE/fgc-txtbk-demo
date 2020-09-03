t = (1:step)*h;
figure
for i = 1:n
plot( t,X(:,3*i-1) - X(:,3*i-2),'LineWidth',2 );
hold on
end
title('sigema of two missles');
xlabel('t/s');
ylabel('rad');
grid on

