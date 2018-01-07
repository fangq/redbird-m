function [Jscatamp,dDdscatamp]=rbjacscatamp(Jd, dcoeff, wavelen, scatpow)

dDdscatamp=-3*dcoeff.*dcoeff*wavelen^(-scatpow);

Jscatamp=Jd.*dDdscatamp;