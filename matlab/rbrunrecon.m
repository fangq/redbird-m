function [dmu, misfit, blocks]=rbrunrecon(niter,cfg,sd,recon,detphi0,f2rid,f2rweight,reform)
%
% [dmu, misfit]=rbreconstep(cfg,sd,recon,detphi0,f2rid,f2rweight, reform)
%
% Perform a single iteration of a Gauss-Newton reconstruction
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the simulation settings stored as a redbird data structure
%     sd: source detector mapping table
%     recon: reconstruction data structure
%     detphi0: measurement data matrix
%     f2rid: the list of element indices of the reconstruction mesh where each forward mesh
%           node is enclosed
%     f2rweight: the barycentric coordinates of the forward mesh nodes inside the
%           reconstruction mesh elements
%
% output:
%     dmua: the unknown update vector computed at the forward mesh nodes 
%          after solving the Gauss-Newton normal equation
%     misfit: the residual betweet the model and the measurement data
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

for iter=1:niter

    [detphi, phi]=rbrunforward(cfg);   % run forward on recon mesh
    if(isfield(cfg,'omega') && cfg.omega>0)
        [Jmua, Jd]=rbfemmatrix(cfg, sd, phi);    % use mex to build Jacobian, 2x faster
    else
        Jmua=rbfemmatrix(cfg, sd, phi);    % use mex to build Jacobian, 2x faster
        %Jmua=rbjacmuafast(sd, phi, cfg.nvol); % use approximated nodal-adjoint for mua
        %Jmua=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol); % or use native code to build nodal-based Jacobian for mua
    end

    if(isa(Jmua,'containers.Map')) % keys are wavelengths, same for Jd
        if(exist('Jd','var'))
            [Jmua,detphi0,detphi,blocks]=rbmultispectral(Jmua, detphi0, detphi, cfg.param, Jd);
        else
            [Jmua,detphi0,detphi,blocks]=rbmultispectral(Jmua, detphi0, detphi, cfg.param);
        end
    else
        if(exist('Jd','var'))
            Jmua=[Jmua, Jd];
            blocks=struct('mua',size(Jmua,2),'dcoeff',size(Jd,2)};
        else
            blocks=struct('mua',size(Jmua,2));
        end

    end

    % from this point on, Jmua is a matrix of the combined Jacobian, unknown
    % blocks are in 

    if(nargin>=8)
        [Jmua,misfit]=rbmatreform(Jmua, detphi0(:), detphi(:), reform);
    else
        misfit=detphi(:)-detphi0(:);
    end
    if(isfield(recon,'elem') && isfield(recon,'node') && nargin>6) % dual-mesh reconstruction
        Jmua_recon=meshremap(Jmua.',f2rid, f2rweight,recon.elem,size(recon.node,1)).'; 
    else % single mesh reconstruction
        Jmua_recon=Jmua;
    end
    lambda=0.1;
    if(isfield(recon,'lambda'))
        lambda=recon.lambda;
    end
    dmu_recon=rbreginv(Jmua_recon, misfit, lambda);  % solve the update on the recon mesh
    if(nargin>6)
        dmu=meshinterp(dmu_recon,f2rid, f2rweight,recon.elem); % interpolate the update to the forward mesh
    else
        dmu=dmu_recon;
    end
end