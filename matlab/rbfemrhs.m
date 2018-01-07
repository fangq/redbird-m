function [rhs,loc,bary,optode]=rbfemrhs(cfg)

% create the right-hand-sides for the FEM system equation, here we solve
% forward systems for both source and detector locations in order to use
% the adjoint method to create Jacobians.

optode=rbgetoptodes(cfg);

if(size(optode,1)<1)
    error('you must provide at least one source or detector');
end

[loc, bary]=tsearchn(cfg.node,cfg.elem,optode);

rhs=zeros(size(cfg.node,1),size(optode,1));
for i=1:size(optode,1)
    if(~isnan(loc(i)))
        rhs(cfg.elem(loc(i),:),i)=bary(i,:)';
    end
end