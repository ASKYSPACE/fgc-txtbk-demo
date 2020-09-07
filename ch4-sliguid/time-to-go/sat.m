function out = sat( s )
%   Saturation function is used to replace sign function to weaken chattering.
%  function inputs  %
%   d    boundary layer thickness
d = 0.1;
if abs(s)<=d
    out = s/d;
else
    out = sign(s);
end
end

