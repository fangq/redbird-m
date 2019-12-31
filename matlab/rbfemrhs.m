function [rhs,loc,bary,optode]=rbfemrhs(cfg)

% create the right-hand-sides for the FEM system equation, here we solve
% forward systems for both source and detector locations in order to use
% the adjoint method to create Jacobians.

optode=rbgetoptodes(cfg);


if(isfield(cfg,'srcpos') && (size(cfg.srcpos,2) == size(cfg.face,1)))
    loc=[];
    bary=[];
    Reff=cfg.reff;
    maxbcnode=max(cfg.face(:));
    srcsum=sum(cfg.srcpos,2);
    srcsum(srcsum==0)=1;

    Adiagbc=cfg.area(:)*((1-Reff)/(9*(1+Reff)));
    Adiagbc=repmat(Adiagbc,1,size(cfg.srcpos,1)).*(cfg.srcpos');

    rhs=sparse(size(cfg.node,1),size(cfg.srcpos,1));
    for i=1:size(cfg.srcpos,1)
        rhs(1:maxbcnode,i)=sparse(cfg.face(:), 1, repmat(Adiagbc(:,i),1,3)/srcsum(i));
    end
    return;
end

if(size(optode,1)<1)
    error('you must provide at least one source or detector');
end

[loc, bary]=tsearchn(cfg.node,cfg.elem,optode);

rhs=sparse(size(cfg.node,1),size(optode,1));
for i=1:size(optode,1)
    if(~isnan(loc(i)))
        rhs(cfg.elem(loc(i),:),i)=bary(i,:)';
    end
end