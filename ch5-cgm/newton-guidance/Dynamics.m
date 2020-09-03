%a             加速度
%v             速度
%the           路径角
%x             横坐标
%y             纵坐标
%X             系统的状态
%U             控制量
%%
%X是当前模型的状态变量，U是模型的控制变量
%返回值为其对应的微分值
%%
function dX = Dynamics(  X ,U)
global m S Cd g rou 
  a = U;  
  v = X(1);
  the = X(2);
  dv = - rou*v^2*S*Cd/2/m-g*sin(the) ; 
  dthe = (a-g*cos(the))/v ;
  dx = v*cos(the) ;  
  dy = v*sin(the) ;
  dX = [dv;dthe;dx;dy] ; 
end