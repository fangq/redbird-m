function [rhs,loc,bary,optode]=rbfemrhs(cfg,sd,wv,md)
%
% newcfg=rbmeshprep(cfg)
%
% Create the right-hand-sides for the FEM system equation, here we solve
% forward systems for both source and detector locations in order to use
% the adjoint method to create Jacobians.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the initial simulation data structure
%
% output:
%     rhs: the right-hand-side of the FEM equation; if multiple sources are
%          used, the rhs is a matrix of dimension Nn x (Ns+Nd), where Nn is 
%          the total number of nodes, and Ns the number of sources and
%          Nd is the total number of detectors
%     loc: the indices of the forward mesh element that encloses each
%          source or detector; Nan means the source is outside of the mesh
%          or a wide-field source/detector
%     bary: the barycentric coordinates of the source/detector if it is
%          enclosed by a tetrahedral element; Nan if outside of the mesh
%     optode: the full source+detector position list returned by
%          rbgetoptodes.m
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if (nargin > 2)
    sd = sd(wv);
    if (nargin > 3 && size(sd,2) > 3)
        sd = sd(sd(:,4) == md,:);
    end
end


[optode,widesrc]=rbgetoptodes(cfg);

if exist('sd','var')
    activeOpt = unique(sd(:,1:2));
else
    activeOpt = size(optode,1) + size(widesrce,1);
end

if(size(optode,1)<1 && size(widesrc,1)<1)
    error('you must provide at least one source or detector');
end

% loc=[];
% bary=[];

if(~isempty(widesrc) && (size(widesrc,2) == size(cfg.node,1)))
    rhs=widesrc.';
    loc=nan*ones(1,size(widesrc,1));
    bary=nan*ones(size(widesrc,1),4);
else
    loc = nan*ones(size(optode,1),1);
    bary = nan*ones(size(optode,1),4);
end

if(isempty(optode))
    return;
end

rhs=sparse(size(cfg.node,1),size(widesrc,1)+size(optode,1));
[newloc, newbary]=tsearchn(cfg.node,cfg.elem,optode(activeOpt,:));

% loc=[loc; newloc];
% bary=[bary; newbary];
loc(activeOpt,:) = newloc;
bary(activeOpt,:) = newbary;

for i=1:size(activeOpt,1)
    if(~isnan(newloc(i)))
%         rhs(cfg.elem(newloc(i),:),i+size(widesrc,1))=newbary(i,:)';
        rhs(cfg.elem(newloc(i),:),activeOpt(i))=newbary(i,:)';
    end
end