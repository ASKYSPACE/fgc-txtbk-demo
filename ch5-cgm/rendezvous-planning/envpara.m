function [df,d2f,r]=envpara(aef,time,emiu)
miu = emiu;  
a  = aef(1);   
e  = aef(2);   
f0 = aef(3);  
n = sqrt(miu/a^3);
E0 = 2*atan2(tan(f0/2),sqrt((1+e)/(1-e)));
tau = (E0 - e*sin(E0))/n;
M = n.*(tau + time);
f = M2f(M,e);
p = a*(1-e^2);
h = sqrt(p*miu);
r = p./(1+e.*cos(f));
df = h./r.^2;
vr  = sqrt(miu/p)*e.*sin(f);
d2f = -2.*vr.*df./r;
function f = M2f(M,e)
M = mod(M,2*pi);
E0 = M;
E = 0;
error = 1e-15;
i = 1;
while (1)
    E = E0 - ((E0-e.*sin(E0))-M)./(1-e.*cos(E0));
    if abs(E-E0)<=error
        break
    else
        E0 = E;
        i = i + 1;
    end
    
end

 f = 2.*atan2(sqrt((1+e)/(1-e)),1./tan(E/2));
