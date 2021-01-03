function [Amat,deldotdel]=rbfemlhs(cfg, deldotdel, wavelength)
%
% [Amat,deldotdel]=rbfemlhs(cfg, deldotdel, wavelength)
%
% create the FEM stiffness matrix (left-hand-side) for solving the diffusion equation
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the initial simulation data structure
%     deldotdel (optional): precomputed operator on the mesh (del_phi dot del_phi)
%         where del represents the gradient; see help rbdeldotdel
%     wavelength (optional): a string or number denoting the wavelength
%
% output:
%     Amat: the left-hand-side matrix of the FEM equation - a sparse matrix
%          of dimension Nn x Nn, where Nn is the number of nodes of the
%          forward mesh
%     deldotdel: if the 2nd input is not given, this function compute
%          deldotdel and return as the 2nd output
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%


nn=size(cfg.node,1);
ne=size(cfg.elem,1);

R_C0=(1./299792458000.);

if(isfield(cfg,'param') && isstruct(cfg.param) && all(structfun(@isempty,cfg.param)==0))
    cfg.prop=rbupdateprop(cfg);
end

prop=cfg.prop;

cfgreff=cfg.reff;
if(isfield(cfg,'mua') && ~isempty(cfg.mua))
    cfgmua=cfg.mua;
end
if(isfield(cfg,'dcoeff') && ~isempty(cfg.dcoeff))
    cfgdcoeff=cfg.dcoeff;
end

if(isa(cfg.prop,'containers.Map'))
    if(nargin<3)
        error('you must specify wavelength');
    end
    if(~ischar(wavelength))
       wavelength=sprintf('%g',wavelength);
    end
    prop=cfg.prop(wavelength);
    cfgreff=cfg.reff(wavelength);
    if(isfield(cfg,'mua') && ~isempty(cfg.mua))
        cfgmua=cfg.mua(wavelength);
    end
    if(isfield(cfg,'dcoeff') && ~isempty(cfg.dcoeff))
        cfgdcoeff=cfg.dcoeff(wavelength);
    end
end

if(nargin>=2)
    if(~isfield(cfg,'mua') || isempty(cfg.mua))
        mua=prop(cfg.seg+1,1);
    else
        mua=cfgmua;
    end
    if(~isfield(cfg,'dcoeff') || isempty(cfg.dcoeff))
        musp=prop(cfg.seg+1,2).*(1-prop(cfg.seg+1,3));
        dcoeff=1./(3*(mua+musp));
    else
        dcoeff=cfgdcoeff;
    end
    if(~isfield(cfg,'nref') || isempty(cfg.nref))
        nref=prop(cfg.seg+1,4);
    else
        nref=cfgreff;
    end
    Reff=cfgreff;

    edges=sort(meshedge(cfg.elem),2);
    
    if(length(mua)==size(cfg.elem,1))
        Aoffd=deldotdel(:,[2:4,6:7,9]).*repmat(dcoeff(:),1,6) + repmat(0.05*mua(:).*cfg.evol(:),1,6);
        Adiag=deldotdel(:,[1,5,8, 10]).*repmat(dcoeff(:),1,4) + repmat(0.10*mua(:).*cfg.evol(:),1,4);
        if(cfg.omega>0)
            Aoffd=complex(Aoffd,repmat(0.05*cfg.omega*R_C0*nref(:).*cfg.evol(:),1,6));
            Adiag=complex(Adiag,repmat(0.10*cfg.omega*R_C0*nref(:).*cfg.evol(:),1,4));
        end
    else
        w1=(1/120)*[2 2 1 1;2 1 2 1; 2 1 1 2;1 2 2 1; 1 2 1 2; 1 1 2 2]';
        w2=(1/60)*(diag([2 2 2 2])+1);
        mua_e=reshape(mua(cfg.elem),size(cfg.elem));
        nref_e=reshape(nref(cfg.elem),size(cfg.elem));
        dcoeff_e=mean(reshape(dcoeff(cfg.elem),size(cfg.elem)),2);
        Aoffd=deldotdel(:,[2:4,6:7,9]).*repmat(dcoeff_e,1,6) + (mua_e*w1).*repmat(cfg.evol(:),1,6);
        Adiag=deldotdel(:,[1,5,8, 10]).*repmat(dcoeff_e,1,4) + (mua_e*w2).*repmat(cfg.evol(:),1,4);
        if(cfg.omega>0)
            Aoffd=complex(Aoffd,(cfg.omega*R_C0)*(nref_e*w1).*repmat(cfg.evol(:),1,6));
            Adiag=complex(Adiag,(cfg.omega*R_C0)*(nref_e*w2).*repmat(cfg.evol(:),1,4));
        end
    end
    edgebc=sort(meshedge(cfg.face),2);
    Adiagbc=cfg.area(:)*((1-Reff)/(12*(1+Reff)));
    Adiagbc=repmat(Adiagbc,1,3);
    Aoffdbc=Adiagbc*0.5;
    
    Amat=sparse([edges(:,1); edges(:,2); cfg.elem(:); edgebc(:,1); edgebc(:,2); cfg.face(:)], ...
                [edges(:,2); edges(:,1); cfg.elem(:); edgebc(:,2); edgebc(:,1); cfg.face(:)], ...
                [Aoffd(:); Aoffd(:); Adiag(:); Aoffdbc(:); Aoffdbc(:); Adiagbc(:)]);
else
    if(size(cfg.elem,2)>4)
        cfg.elem(:,5:end)=[];
    end
    [Adiag, Aoffd, deldotdel]=rbfemmatrix(cfg);
    Amat = sparse([cfg.rows,cfg.cols,(1:nn)],[cfg.cols,cfg.rows,(1:nn)],[Aoffd,Aoffd,Adiag],nn,nn);
    deldotdel=deldotdel';
end
