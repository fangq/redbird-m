function [dmu, misfit]=rbreconstep(cfg,sd,recon,detphi0,f2rid,f2rweight,reform)
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

[detphi, phi]=rbrunforward(cfg);   % run forward on recon mesh
%Jmua=rbfemmatrix(cfg, sd, phi);    % use mex to build Jacobian, 2x faster
Jmua=rbjacmuafast(sd, phi, cfg.nvol); % use approximated nodal-adjoint for mua
%Jmua=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol); % or use native code to build nodal-based Jacobian for mua

if(isa(Jmua,'containers.Map'))
    [Jmua,detphi0,detphi]=rbmultispectral(Jmua, detphi0, detphi);
end
if(nargin>=7)
    [Jmua,misfit]=rbmatreform(Jmua, detphi0(:), detphi(:),reform);
else
    misfit=detphi(:)-detphi0(:);
end
if(isfield(recon,'elem') && isfield(recon,'node') && nargin>5) % dual-mesh reconstruction
    Jmua_recon=meshremap(Jmua.',f2rid, f2rweight,recon.elem,size(recon.node,1)).'; 
else % single mesh reconstruction
    Jmua_recon=Jmua;
end
lambda=0.1;
if(isfield(recon,'lambda'))
    lambda=recon.lambda;
end
dmu_recon=rbreginv(Jmua_recon, misfit, lambda);  % solve the update on the recon mesh
if(nargin>5)
    dmu=meshinterp(dmu_recon,f2rid, f2rweight,recon.elem); % interpolate the update to the forward mesh
else
    dmu=dmu_recon;
end
