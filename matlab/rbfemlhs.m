function [Amat,deldotdel]=rbfemlhs(cfg, deldotdel)

% create the FEM stiffness matrix (left-hand-side) for solving the diffusion equation

nn=size(cfg.node,1);
ne=size(cfg.elem,1);

if(nargin==2)
    if(~isfield(cfg,'mua') || isempty(cfg.mua))
        mua=cfg.prop(cfg.elemprop+1,1);
    else
        mua=cfg.mua;
        if(length(mua)==nn)
            mua=mean(reshape(mua(cfg.elem),ne,size(cfg.elem,2)),2);
        end
    end
    if(~isfield(cfg,'dcoeff') || isempty(cfg.dcoeff))
        musp=cfg.prop(cfg.elemprop+1,2).*(1-cfg.prop(cfg.elemprop+1,3));
        dcoeff=1./(3*(mua+musp));
    else
        dcoeff=cfg.dcoeff;
        if(length(dcoeff)==nn)
            dcoeff=mean(reshape(dcoeff(cfg.elem),ne,size(cfg.elem,2)),2);
        end
    end
    if(~isfield(cfg,'nref') || isempty(cfg.nref))
        nref=cfg.prop(cfg.elemprop+1,4);
    else
        nref=cfg.nref;
        if(length(nref)==nn)
            nref=mean(reshape(nref(cfg.elem),ne,size(cfg.elem,2)),2);
        end
    end
    Reff=cfg.reff;

    edges=sort(meshedge(cfg.elem),2);
    Aoffd=deldotdel(:,[2:4,6:7,9]).*repmat(dcoeff(:),1,6) + repmat(0.05*mua(:).*cfg.evol(:),1,6);
    Adiag=deldotdel(:,[1,5,8, 10]).*repmat(dcoeff(:),1,4) + repmat(0.10*mua(:).*cfg.evol(:),1,4);
    if(cfg.omega>0)
        R_C0=(1./299792458000.);
        Aoffd=complex(Aoffd,repmat(0.05*cfg.omega*R_C0*nref(:).*cfg.evol(:),1,6));
        Adiag=complex(Adiag,repmat(0.10*cfg.omega*R_C0*nref(:).*cfg.evol(:),1,4));
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
