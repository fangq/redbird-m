function [Jscatamp,dDdscatpow]=rbjacscatpow(Jd, dcoeff, wavelen)

dDdscatpow=dcoeff*log(wavelen);

Jscatamp=Jd.*dDdscatpow;