function Jchrome=rbjacchrome(Jmua, chromorphores)
%
% Jchrome=rbjacchrome(Jmua, extin)
%
% Building the Jacobian matrices for chromorphores from Jacobian of the 
% absorption coeff and extinction coeff.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jmua: the Jacobian for absorption mua, must be a containers.Map
%     chromorphores: list of chmorphores supported by rbextinction
%
% output:
%     Jchrome: the Jacobian for all chmorphores, in a struct Jchrome.{hbo,hbr,...} 
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(~ismatrix(Jmua)) % if a containers.Map, flatten it
    error('Jmua must be a containers.Map with elements of each wavelength');
end

wavelengths=keys(Jmua);
extin=rbextinction(wavelengths, chromorphores);

Jchrome=struct;
for i=1:length(chromorphores)
    Jchrome.(chromorphores{i})=rbmatflat(Jmua,extin(:,i));
end

% J=[J(chrom1),J(chrom2),...,J(chromM)]