function [dmu, misfit]=rbreconstep(cfg,sd,recon,f2rid,f2rweight)

[detphi, phi]=rbrunforward(cfg);   % run forward on recon mesh
Jmua=rbfemmatrix(cfg, sd, phi);    % use mex to build Jacobian, 2x faster
%Jmua=rbjacmuafast(sd, phi, cfg.nvol); % use approximated nodal-adjoint for mua
%Jmua=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol); % or use native code to build nodal-based Jacobian for mua
Jmua_recon=meshremap(Jmua.',f2rid, f2rweight,recon.elem,size(recon.node,1)).'; 
misfit=detphi0(:)-detphi(:);       % calculate data-misfit
[Jmua_recon,misfit]=rbmatreform(Jmua_recon, detphi0(:), detphi(:), 'logphase');
resid(i)=sum(abs(misfit));         % store the residual
dmu_recon=rbreginv(Jmua_recon, misfit, 0.05);  % solve the update on the recon mesh
dmu=meshinterp(dmu_recon,f2rid, f2rweight,recon.elem); % interpolate the update to the forward mesh