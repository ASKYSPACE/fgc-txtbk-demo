function obsp = obsevol(obs_s0,time,aef,emiu,opts,flag)

if flag == 0
obsp = kron([obs_s0(1) obs_s0(2) obs_s0(3)],ones(length(time),1));
elseif flag == 1
    span = linspace(0,time(end),1E4);
    [t,obsp_t] = ode45(@FuncT_H,span,obs_s0,opts,aef,emiu);
    obsp = interp1(t,obsp_t,time,'spline');
else
    error('flag=1 or flag =0');
end
function ds=FuncT_H(t,s,aef,emiu)
ds = zeros(6,1);
x  = s(1);
y  = s(2);
z  = s(3);
vx = s(4);
vy = s(5);
vz = s(6);
[dsit,ddsit,r] = envpara(aef,t,emiu);
miu = emiu;
ds(1) = vx;
ds(2) = vy;
ds(3) = vz;
ds(4) = (dsit^2 + 2*miu/(r^3))*x + ddsit*y + 2*dsit*vy;
ds(5) = -ddsit*x+(dsit^2-miu/(r^3))*y-2*dsit*vx;
ds(6) = -miu*z/(r^3);
