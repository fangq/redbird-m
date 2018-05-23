function [Amat,deldotdel]=rbfemlhs(cfg, deldotdel)

% create the FEM stiffness matrix (left-hand-side) for solving the diffusion equation

nn=size(cfg.node,1);

if(nargin==2)
    mua=cfg.prop(cfg.elemprop+1,1);
    musp=cfg.prop(cfg.elemprop+1,2).*(1-cfg.prop(cfg.elemprop+1,3));
    dcoeff=1./(3*(mua+musp));
    nref=cfg.prop(cfg.elemprop+1,4);

    edges=sort(meshedge(cfg.elem),2);
    Aoffd=deldotdel(:,[2:4,6:7,9]).*repmat(dcoeff(:),1,6) + repmat(0.05*mua(:).*cfg.evol(:),1,6);
    Adiag=deldotdel(:,[1,5,8, 10]).*repmat(dcoeff(:),1,4) + repmat(0.10*mua(:).*cfg.evol(:),1,4);
    if(cfg.omega>0)
        R_C0=(1./299792458000.);
        Aoffd=complex(Aoffd,repmat(0.05*cfg.omega*R_C0*nref(:).*cfg.evol(:),1,6));
        Adiag=complex(Adiag,repmat(0.10*cfg.omega*R_C0*nref(:).*cfg.evol(:),1,4));
    end

    edgebc=sort(meshedge(cfg.face),2);
    Reff=rbgetreff(cfg.prop(cfg.elemprop(1)+1,4), cfg.prop(1,4));
    Adiagbc=cfg.area(:)*(1-Reff)/(12*(1+Reff));
    Adiagbc=repmat(Adiagbc,1,3);
    Aoffdbc=Adiagbc*0.5;
    
    Amat=sparse([edges(:,1); edges(:,2); cfg.elem(:); edgebc(:,1); edgebc(:,2); cfg.face(:)], ...
                [edges(:,2); edges(:,1); cfg.elem(:); edgebc(:,2); edgebc(:,1); cfg.face(:)], ...
                [Aoffd(:); Aoffd(:); Adiag(:); Aoffdbc(:); Aoffdbc(:); Adiagbc(:)]);
else
    [Adiag, Aoffd, deldotdel]=rbfemmatrix(cfg);
    Amat = sparse([cfg.rows,cfg.cols,(1:nn)],[cfg.cols,cfg.rows,(1:nn)],[Aoffd,Aoffd,Adiag],nn,nn);
    deldotdel=deldotdel';
end
