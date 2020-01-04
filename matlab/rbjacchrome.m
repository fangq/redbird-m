function Jchrome=rbjacchrome(Jmua, extin)
%
% Jchrome=rbjacchrome(Jmua, extin)
%
% Building the Jacobian matrices for chromorphores from Jacobian of the 
% absorption coeff and extinction coeff.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jmua: the Jacobian for absorption coeff. mua
%     extin: the extinction coeff as a matrix: rows: wavelengths; columns: chromorphores
%
% output:
%     Jmua_n: the nodal Jacobian for absorption coeff. mua
%     Jmua_e: the element-wise Jacobian for absorption coeff. mua
%     Jd_n: (optional) the nodal Jacobian for diffusion coeff D
%     Jd_e: (optional) the element-wise Jacobian for diffusion coeff D
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

extin=extin';
extin=extin(:)';
if(~ismatrix(Jmua))
    Jmua=rbmatflat(Jmua); %J=[J(wavelen1); J(wavelen2) ;...; J(wavelenN)]
end
Jchrome=kron(extin, Jmua);

% J=[J(chrom1),J(chrom2),...,J(chromM)]