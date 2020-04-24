function [newJ, newy0, newphi]=rbinvmatrix(Jacob, y0, phi)
%
% [newJ, newy0, newphi]=rbinvmatrix(Jacob, y0, phi)
%
% Concatenate multi-special forward modeling data into a single linear
% system
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jacob: the multispectral Jacobian as a containers.Map object-each
%         wavelength is a key
%     y0: the multispectral measurement data 
%     phi: the forward solution at multiple wavelengths
%
% output:
%     newJ: the contatenated Jacobian - dimension Nw*Ns*Nd x Nn, where
%          Nw - the number of wavelengths
%          Ns - the number of sources
%          Nd - the number of detectors
%          Nn - the number of nodes
%     newy0: the contatenated measurement data vector
%     newphi: the contatenated forward simulation of measurements
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

newJ=[];
newy0=[];
newphi=[];

if(isa(Jacob,'containers.Map') && isa(y0,'containers.Map') && isa(phi,'containers.Map'))
    wv=keys(Jacob);
    if(~all(ismember(wv,keys(y0))) || ~all(ismember(keys(phi),keys(y0))))
        error('Jacob, y0 and phi must share the same key');
    end
    for i=1:length(wv)
        newJ=[newJ; Jacob(wv{i})];
        newy0=[newy0; y0(wv{i})];
        newphi=[newphi; phi(wv{i})];
    end
end