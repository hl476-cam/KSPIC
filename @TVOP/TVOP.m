function  res = TVOP(minus1D,TVtype,aTV)

%res = TVOP()
%
% Implements a spatial finite-differencing operator.
%
% (c) Michael Lustig 2007

res.adjoint = 0;
res.minus1D = 0;
res.TVtype = 0;
res.weight = 0;

if(nargin)
    res.minus1D = minus1D;
    if(exist('TVtype','var'))
        res.TVtype = TVtype;
    end
    if(exist('aTV','var'))
        res.weight = aTV;
    end
end 

res = class(res,'TVOP');

