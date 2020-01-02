function [rhs,loc,bary,optode]=rbfemrhs(cfg)

% create the right-hand-sides for the FEM system equation, here we solve
% forward systems for both source and detector locations in order to use
% the adjoint method to create Jacobians.

[optode,widesrc]=rbgetoptodes(cfg);

if(size(optode,1)<1 && size(widesrc,1)<1)
    error('you must provide at least one source or detector');
end

loc=[];
bary=[];
rhs=sparse(size(cfg.node,1),size(widesrc,1)+size(optode,1));

if(~isempty(widesrc) && (size(widesrc,2) == size(cfg.face,1)))
    Reff=cfg.reff;
    maxbcnode=max(cfg.face(:));

    Adiagbc=cfg.area(:)*((1-Reff)/(18*(1+Reff)));
    Adiagbc=repmat(Adiagbc,1,size(widesrc,1)).*(widesrc');

    for i=1:size(widesrc,1)
        rhs(1:maxbcnode,i)=sparse(cfg.face(:), 1, repmat(Adiagbc(:,i),1,3));
        rhs(:,i)=rhs(:,i)/sum(rhs(:,i));
    end
    loc=nan*ones(1,size(widesrc,1));
    bary=nan*ones(size(widesrc,1),4);
end

if(isempty(optode))
    return;
end

[newloc, newbary]=tsearchn(cfg.node,cfg.elem,optode);

loc=[loc; newloc];
bary=[bary; newbary];

for i=1:size(optode,1)
    if(~isnan(newloc(i)))
        rhs(cfg.elem(newloc(i),:),i+size(widesrc,1))=newbary(i,:)';
    end
end