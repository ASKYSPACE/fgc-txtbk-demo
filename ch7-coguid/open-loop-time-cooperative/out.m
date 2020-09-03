global n step
P = zeros(10,n,2);
for i = 1:step
    for j = 1:n
    P(i,j,:) =Pt - X( i,3*j )*[ cos( X( i,3*j-2 ) ),sin( X( i,3*j-2 ) ) ];
    end
end
figure
for j = 1:n
plot( P(:,j,1),P(:,j,2),'LineWidth',2 );
grid on
hold on
end
title('µ¼µ¯¹ì¼£');
xlabel('XÖá/Ã×');
ylabel('YÖá/Ã×');