function [dmu, misfit]=rbreconstep(cfg,sd,recon,detphi0,f2rid,f2rweight, reform)

[detphi, phi]=rbrunforward(cfg);   % run forward on recon mesh
%Jmua=rbfemmatrix(cfg, sd, phi);    % use mex to build Jacobian, 2x faster
Jmua=rbjacmuafast(sd, phi, cfg.nvol); % use approximated nodal-adjoint for mua
%Jmua=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol); % or use native code to build nodal-based Jacobian for mua
[Jmua,detphi0,detphi]=rbinvmatrix(Jmua, detphi0, detphi);
if(nargin>=7)
    [Jmua,misfit]=rbmatreform(Jmua, detphi0(:), detphi(:),reform);
end
Jmua_recon=meshremap(Jmua.',f2rid, f2rweight,recon.elem,size(recon.node,1)).'; 
dmu_recon=rbreginv(Jmua_recon, misfit, 0.05);  % solve the update on the recon mesh
dmu=meshinterp(dmu_recon,f2rid, f2rweight,recon.elem); % interpolate the update to the forward mesh