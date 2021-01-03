function varargout=rbrun(cfg,recon,varargin)
%
% phi=rbrecon(cfg)
% prop=rbrecon(cfg,recon)
%   or
% prop=rbrecon(cfg,recon,mode)
%
% Perform reconstruction by fitting optical properties with the data vector
% provided recon or recon.data
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird data structure, if cfg is the only input, rbrun is
%            the same as rbrunforward
%     recon: if recon is a vector, it must be the measurement data for all
%               src/detector pairs
%            if recon is a struct, it can have the following subfields
%               recon.data (required): the data to be fitted for
%               recon.sdmap: the sd map of the data vector
%               recon.node, recon.elem: the reconstruction mesh
%               recon.iter: number of iterations, if not present, use 10
%               recon.lambda: Tikhonov regularization parameter, if not
%                     present, use 0.1
%               recon.prop0 or recon.param0: initial guess of mua/mus or
%                     parameters
%     mode: if mode='image' or 'r' (default), rbrecon fit the data in recon
%           if mode='simu' or 'x', rbrecon first creates simulated forward
%               data and then reconstruct using the recon settings
%           if mode='bulk' or 'b', bulk fitting of a single segment
%                    if recon.seg presents, this mode returns optical
%                    properties per segment
%           if mode='fuzzy' or 'f', bulk fitting of a labeled segmentation
%                    recon.seg must present, this mode returns optical
%                    properties per segment by backmapping optical
%                    properties back to the volume
%           if mode='prior' or 'p', compositional-prior reconstruction
%           if mode='roiprior', compositional-prior reconstruction with
%                    roi and background probabalistic maps
%           if mode='roionly', 2-compositional-prior reconstruction with
%                    roi/non-roi probabalistic map
%
% output:
%     prop: reconstructed optical properties or parameters
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

mode='image';
if(nargin>2)
    mode=varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run forward for the heterogeneous domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin==1)
    varargout{1:nargout}=rbrunforward(cfg);
    return;
elseif(isstruct(recon) && isfield(recon,'srcpos'))
    detphi0=rbrunforward(recon);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Reset the domain to a homogeneous medium for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg=rbsetmesh(cfg,recon.node,recon.elem,recon.prop,ones(size(node,1),1));

if(isfield(recon,'sdmap'))
    sd=recon.sdmap;
else
    sd=rbsdmap(cfg);
end

[f2rid, f2rweight]=tsearchn(recon.node,recon.elem,cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run 10 iterations to recover mua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxiter=10;
if(isfield(recon,'iter'))
    maxiter=iter;
end
resid=zeros(1,maxiter);

cfg.mua=ones(size(cfg.node,1),1)*cfg.prop(cfg.seg(1)+1,1);

for i=1:maxiter
    tic
    [detphi, phi]=rbrunforward(cfg);   % run forward on recon mesh
    Jmua=rbfemmatrix(cfg, sd, phi);    % use mex to build Jacobian, 2x faster
    %Jmua=rbjacmuafast(sd, phi, cfg.nvol); % use approximated nodal-adjoint for mua
    %Jmua=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol); % or use native code to build nodal-based Jacobian for mua
    Jmua_recon=meshremap(Jmua.',f2rid, f2rweight,recon.elem,size(recon.node,1)).'; 
    misfit=detphi0(:)-detphi(:);       % calculate data-misfit
    [Jmua_recon,misfit]=rbcreateinv(Jmua_recon, detphi0(:), detphi(:), 'logphase');
    resid(i)=sum(abs(misfit));         % store the residual
    dmu_recon=rbreginv(Jmua_recon, misfit, 0.05);  % solve the update on the recon mesh
    dmu=meshinterp(dmu_recon,f2rid, f2rweight,recon.elem); % interpolate the update to the forward mesh
    cfg.mua=cfg.mua + dmu(:);          % update forward mesh mua vector
    fprintf(1,'iter [%4d]: residual=%e, relres=%e (time=%f s)\n',i, resid(i), resid(i)/resid(1), toc);
end
