function [Jscatamp,dDdscatamp]=rbjacscatamp(Jd, dcoeff, wavelen, scatpow)
%
% [Jscatamp,dDdscatamp]=rbjacscatamp(Jd, dcoeff, wavelen, scatpow)
%
% Create the Jacobian matrix for the scattering amplitude using the Jacobian 
% of the diffusion coeff (D)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jd: the jacobian of the diffusion coeff.
%     dcoeff: the diffusion coefficient values at each node
%     wavelen: the list of wavelengths
%     scatpow: the scattering power of the current estimate of scattering power
%
% output:
%     Jscatamp: the Jacobian of the scattering amplitude parameter
%     dDdscatamp: partial derivative of D - diffusion coeff - to the scat amplitude
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

dDdscatamp=-3*dcoeff.*dcoeff*(wavelen/500)^(-scatpow);

Jscatamp=Jd.*dDdscatamp;