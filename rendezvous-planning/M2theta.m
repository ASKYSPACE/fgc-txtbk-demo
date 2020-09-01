function f = M2theta(M,e)
E0 = M;
E = 0;
error = 1e-3;
i = 1;
while (1)
    E = E0 - ((E0-e.*sin(E0))-M)./(1-e.*cos(E0));
    if abs(E-E0)<=error
        break
    else
        E0 = E;i = i+1
    end
end
 f = 2.*atan2(sqrt((1+e)/(1-e)),1./tan(E/2));