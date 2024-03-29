function [newJ, newy0, newphi]=rbmultispectral(Jmua, y0, phi, params, Jd, prop)
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
%     prop: must be a struct as cfg.prop
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
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

newJ=struct;
newy0=[];
newphi=[];

if(isa(Jmua,'containers.Map'))
    wv=keys(Jmua);
    paramlist=fieldnames(params);
    if(nargin>5 && length(intersect(paramlist,{'scatamp','scatpow'}))==2)
        dcoeff = containers.Map();
        for i=1:length(wv)
            dtemp=prop(wv);
            if(size(dtemp,1)<size(Jd(wv{i}),2)) % label based
                dtemp=dtemp(2:end,:);
            end
            dcoeff(wv{i})=1/(3*(dtemp(:,1)+dtemp(:,2)));
        end
        Jscat=rbjacscat(Jd, dcoeff, params.scatpow, wv);
    end
    chromophores=intersect(paramlist,{'hbo','hbr','water','lipids','aa3'});

    newJ=rbjacchrome(Jmua,chromophores);
    if(exist('Jscat','var') && isstruct(Jscat))
        allkeys=fieldnames(Jscat);
        for i=1:length(allkeys)
            newJ.(allkeys{i})=Jscat.(allkeys{i});
        end
        clear Jscat;
    end
else
    newJ.mua=Jmua;
    if(nargin>4 && ~isa(Jd,'containers.Map'))
        newJ.dcoeff=Jd;
    end
end

if(nargin>4 && ~isa(Jd,'containers.Map'))
    newJ.dcoeff=Jd;
end
if(~isa(y0,'containers.Map'))
    newy0=y0;
else
    wv=keys(y0);
    for i=1:length(wv)
        newy0=[newy0; reshape(y0(wv{i}),[],1)];
    end
end
if(~isa(phi,'containers.Map'))
    newphi=phi;
else
    wv=keys(phi);
    for i=1:length(wv)
        newphi=[newphi; reshape(phi(wv{i}),[],1)];
    end
end
