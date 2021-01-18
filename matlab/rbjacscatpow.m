function [Jscatpow,dDdscatpow]=rbjacscatpow(Jd, dcoeff, wavelen)
%
% [Jscatpow,dDdscatpow]=rbjacscatpow(Jd, dcoeff, wavelen)
%
% Create the Jacobian matrix for the scattering power using the Jacobian 
% of the diffusion coeff (D)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jd: the jacobian of the diffusion coeff.
%     dcoeff: the diffusion coefficient values at each node
%     wavelen: the list of wavelengths
%
% output:
%     Jscatpow: the Jacobian of the scattering power parameter
%     dDdscatpow: partial derivative of D - diffusion coeff - to the scat power
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

dDdscatpow=dcoeff*log(wavelen/500);

Jscatpow=Jd.*dDdscatpow;