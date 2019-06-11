function [Jscatpow,dDdscatpow]=rbjacscatpow(Jd, dcoeff, wavelen)

dDdscatpow=dcoeff*log(wavelen);

Jscatpow=Jd.*dDdscatpow;