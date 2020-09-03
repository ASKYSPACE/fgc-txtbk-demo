for s = 1:12
errorvx(s) = max(abs(vxk(s+1,:)-vxk(s,:)));
errorvy(s) = max(abs(vyk(s+1,:)-vyk(s,:)));
errorvz(s) = max(abs(vzk(s+1,:)-vzk(s,:))); 
end