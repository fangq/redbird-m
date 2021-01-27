function [Amat,deldotdel]=rbfemlhs(cfg, deldotdel, wavelength)
%
% [Amat,deldotdel]=rbfemlhs(cfg)
%   or
% [Amat,deldotdel]=rbfemlhs(cfg, wavelength)
% [Amat,deldotdel]=rbfemlhs(cfg, deldotdel, wavelength)
%
% create the FEM stiffness matrix (left-hand-side) for solving the
% diffusion equation at a given wavelength
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

% cfg.prop is updated from cfg.param and contains the updated mua musp.
% if cfg.param is node/elem based, cfg.prop is updated to have 4 columns
% with mua/musp being the first two columns
% if cfg.param is segmentation based, cfg.prop has the same format as mcx's
% prop, where the first row is label 0, and total length is Nseg+1

prop=cfg.prop;
cfgreff=cfg.reff;
omega=cfg.omega;

if(nargin==2 && numel(deldotdel)==1)
    wavelength=deldotdel;
end

if(isa(cfg.prop,'containers.Map')) % if multiple wavelengths, take current
    if(nargin<3)
        error('you must specify wavelength');
    end
    if(~ischar(wavelength))
       wavelength=sprintf('%g',wavelength);
    end
    prop=cfg.prop(wavelength);
    cfgreff=cfg.reff(wavelength);
    if(isa(omega,'containers.Map'))
       omega=omega(wavelength);
    end
end

% if deldotdel is provided, call native code; otherwise, call mex

if(nargin>=2 && numel(deldotdel)>1)
    % get mua from prop(:,1) if cfg.prop has wavelengths
    if(size(prop,1)==nn || size(prop,1)==ne)
        mua=prop(:,1);
        if(size(prop,2)<3)
            musp=prop(:,2);
        else
            musp=prop(:,2).*(1-prop(:,3));
        end
    elseif(size(prop,1)<min([nn ne])) % use segmentation based prop list
        mua=prop(cfg.seg+1,1);
        if(size(prop,2)<3)
            musp=prop(cfg.seg+1,2); % assume g is 0
        else
            musp=prop(cfg.seg+1,2).*(1-prop(cfg.seg+1,3));
        end
    end
    dcoeff=1./(3*(mua+musp));

    if(isfield(cfg,'bulk') && isfield(cfg.bulk,'n'))
        nref=cfg.bulk.n;
    elseif(isfield(cfg,'seg') && size(prop,1)<min([nn,ne]))
        nref=rbgetbulk(cfg);
        if(isa(nref,'containers.Map'))
            nref=nref(wavelength);
        end
        nref=nref(4);
    else
        nref=prop(:,4);
    end
    Reff=cfgreff;

    edges=sort(meshedge(cfg.elem),2);
    
    % what LHS matrix needs is dcoeff and mua, must be node or elem long
    if(length(mua)==size(cfg.elem,1))  % element based property
        Aoffd=deldotdel(:,[2:4,6:7,9]).*repmat(dcoeff(:),1,6) + repmat(0.05*mua(:).*cfg.evol(:),1,6);
        Adiag=deldotdel(:,[1,5,8, 10]).*repmat(dcoeff(:),1,4) + repmat(0.10*mua(:).*cfg.evol(:),1,4);
        if(omega>0)
            Aoffd=complex(Aoffd,repmat(0.05*omega*R_C0*nref(:).*cfg.evol(:),1,6));
            Adiag=complex(Adiag,repmat(0.10*omega*R_C0*nref(:).*cfg.evol(:),1,4));
        end
    else  % node based properties
        w1=(1/120)*[2 2 1 1;2 1 2 1; 2 1 1 2;1 2 2 1; 1 2 1 2; 1 1 2 2]';
        w2=(1/60)*(diag([2 2 2 2])+1);
        mua_e=reshape(mua(cfg.elem),size(cfg.elem));
        if(length(nref)==1)
            nref_e=nref*ones(size(cfg.elem));
        else
            nref_e=reshape(nref(cfg.elem),size(cfg.elem));
        end
        dcoeff_e=mean(reshape(dcoeff(cfg.elem),size(cfg.elem)),2);
        Aoffd=deldotdel(:,[2:4,6:7,9]).*repmat(dcoeff_e,1,6) + (mua_e*w1).*repmat(cfg.evol(:),1,6);
        Adiag=deldotdel(:,[1,5,8, 10]).*repmat(dcoeff_e,1,4) + (mua_e*w2).*repmat(cfg.evol(:),1,4);
        if(cfg.omega>0)
            Aoffd=complex(Aoffd,(omega*R_C0)*(nref_e*w1).*repmat(cfg.evol(:),1,6));
            Adiag=complex(Adiag,(omega*R_C0)*(nref_e*w2).*repmat(cfg.evol(:),1,4));
        end
    end
    % add partial current boundary condition
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
    cfg.prop=prop; % use property of the current wavelength
    [Adiag, Aoffd, deldotdel]=rbfemmatrix(cfg);
    Amat = sparse([cfg.rows,cfg.cols,(1:nn)],[cfg.cols,cfg.rows,(1:nn)],[Aoffd,Aoffd,Adiag],nn,nn);
    deldotdel=deldotdel';
end
