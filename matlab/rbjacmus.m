function Jmus=rbjacmus(Jd, mus, g)

if(nargin<3)
    g=0;
end
Jmus=-Jd/(3*mus*mus*(1-g));