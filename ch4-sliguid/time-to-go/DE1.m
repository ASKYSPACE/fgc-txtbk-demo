function dy = DE1( ~,y,V_m,c,k,epsi,lamda_d )
%  Sliding mode variable structure guidance law based on general exponential 
%  reaching law is used for trajectory integration.
%  There are 8 integral variables. 
% r           relative distance
% lamda       line of sight angle
% theta_t     flight path angle of target
% theta_m     flight path angle of missile
% x_m         horizontal distance
% H_m         height
% d_lamda     derivative of line of sight angle
% s           sliding face
dy = zeros(8,1);
% Stationary target
a_t = 0;
V_t = 0;
% Differential function
eta_m = y(2)-y(4);
eta_t = y(2)-y(3);
dy(1) = -V_m*cos(eta_m)+V_t*cos(eta_t);
dy(2) = y(7);
s = dy(2)+c*(y(2)-lamda_d);  % sliding face
s = y(8);
a_c = ( y(1)*(-2*dy(1)/y(1)+c+k)*dy(2)+c*y(1)*k*(y(2)-lamda_d)+...
    a_t*cos(eta_t)+epsi*y(1)*sat(s) )/cos(eta_m);
dy(8) = -k*y(8)-epsi*sat(y(8));  % reaching law
a_m = a_c;
dy(7) = -2*dy(1)/y(1)*y(7)-cos(eta_m)/y(1)*a_m+cos(eta_t)/y(1)*a_t;
dy(3) = 0;
dy(4) = a_m/V_m;
dy(5) = V_m*cos(y(4));
dy(6) = V_m*sin(y(4));
end