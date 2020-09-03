function y = RK_4(X)
global h step v n acc acc_b N
% X = X';
K1 = h*solve( X );
K2 = h*solve( X+K1/2);
K3 = h*solve( X+K2/2);
K4 = h*solve( X+K3);
y = X+( K1 + 2*K2 + 2*K3 + K4 )/6;
 X_1 = (y-X)/h;
 for i = 1:n
 acc(step,i) = v(i)*X_1(3*i-1);  %导弹i的加速度
 acc_b(step,i) = acc(step,i) - N*v(i)*X_1(3*i-2);  %导弹i的附加加速度
 end
% n_m(step,2) = v(2)*X_1(5);  %导弹2的加速度
% a(step,1) = v(1)*( X_1(2) - 3*X_1(1) );  % 附加加速度
% a(step,2) = v(2)*( X_1(5) - 3*X_1(4) );  % 附加加速度
% a(step,3) = -v(1)*sin( X(2)-X(1) )/X(3);
% a(step,4) = -v(2)*sin( X(5)-X(4) )/X(6);
% a(step,3) = v(1)*( 3*X_1(1) );  % 原始加速度
% a(step,4) = v(2)*( 3*X_1(4) );  % 原始加速度



