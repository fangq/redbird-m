function [newJ, newy0, newphi, blocks]=rbmultispectral(Jmua, y0, phi, params, Jd)
%
% [newJ, newy0, newphi]=rbmultispectral(Jmua, y0, phi, paramlist, Jd)
%
% Concatenate multi-special forward modeling data into a single linear
% system
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jmua: the multispectral Jacobian as a containers.Map object with each
%         wavelength as a key
%     y0: the multispectral measurement data 
%     phi: the forward solution at multiple wavelengths
%     params: must be a struct as cfg.param
%
% output:
%     newJ: the contatenated Jacobian - dimension (Nw*Ns*Nd) x ((Nn*Np)+Nn*2), where
%          Nw - the number of wavelengths
%          Ns - the number of sources
%          Nd - the number of detectors
%          Nn - the number of nodes
%          Np - the number of parameter species
%          Nn*2 stores the scat-amplitude and scat-power if Jd is given
%     newy0: the contatenated measurement data vector
%     newphi: the contatenated forward simulation of measurements
%     blocks: a struct recording the names and lengths of the unknown
%          blocks
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

newJ=[];
newy0=[];
newphi=[];

if(isa(Jmua,'containers.Map') && isa(y0,'containers.Map') && isa(phi,'containers.Map'))
    wv=keys(Jmua);
    if(~all(ismember(wv,keys(y0))) || ~all(ismember(keys(phi),keys(y0))))
        error('Jacob, y0 and phi must share the same key');
    end
    paramlist=fieldnames(params);
    if(nargin>4 && length(intersect(paramlist,{'scatamp','scatpow'}))==2)
        for i=1:length(wv)
            newJ=rbjacscat(Jd, dcoeff, params.scatpow, wv);
        end
        blocknames=struct('scatamp',size(Jd,2),'scatpow',size(Jd,2));
    end
    chromophores=intersect(paramlist,{'hbo','hbr','water','lipids','aa3'});

    newJ=[rbjacchrome(Jmua,rbextinction(wv, chromophores)) , newJ];
    blocknames=[chromophores,blocknames];
    blocknames=struct('scatamp',size(Jd,2),'scatpow',size(Jd,2));

    for i=1:length(wv)
        newy0=[newy0; reshape(y0(wv{i}),[],1)];
        newphi=[newphi; reshape(phi(wv{i}),[],1)];
    end
else
    if(~isa(Jmua,'containers.Map'))
        newJ=Jmua;
        blocknames=struct('mua',size(Jmua,2));
    end
    if(nargin>4 && ~isa(Jd,'containers.Map'))
        newJ=[Jmua Jd];
        blocknames=struct('mua',size(Jmua,2),'dcoeff',size(Jd,2));
    end
    if(~isa(y0,'containers.Map'))
        newy0=y0;
    end
    if(~isa(phi,'containers.Map'))
        newphi=phi;
    end
end