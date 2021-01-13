%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography, 
%      Copyright Qianqina Fang, 2018
%
% In this example, we show the most basic usage of Redbird.
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(fullfile(pwd, '../matlab'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg cfg0

s0=[70, 50, 20];

[nobbx,fcbbx]=meshabox([40 0 0], [160, 120, 60], 10);
[nosp,fcsp]=meshasphere(s0, 5, 1);
[no,fc]=mergemesh(nobbx, fcbbx, nosp, fcsp);

[cfg0.node, cfg0.elem]=s2m(no,fc(:,1:3),1,40,'tetgen',[41 1 1;s0]);

%[cfg.node, cfg.face, cfg.elem]=meshabox([0 0 0],[60 60 30],3);
nn=size(cfg0.node,1);
cfg0.seg=cfg0.elem(:,5);
cfg0.srcdir=[0 0 1];

[xi,yi]=meshgrid(60:20:140,20:20:100);
cfg0.srcpos=[xi(:),yi(:),zeros(numel(yi),1)];
cfg0.detpos=[xi(:),yi(:),60*ones(numel(yi),1)];

cfg0.prop=[
    0 0 1 1
    0.008 1 0 1.37
    0.016 1 0 1.37
];

z0=1/(cfg0.prop(2,1)+cfg0.prop(2,2)*(1-cfg0.prop(2,3)));

cfg0.srcpos(:,3)=cfg0.srcpos(:,3)+z0;
cfg0.detpos(:,3)=cfg0.detpos(:,3)-z0;

cfg0.omega=2*pi*70e6;
cfg0.omega=0;

cfg=cfg0;

cfg0=rbmeshprep(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run forward for the heterogeneous domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detphi0=rbrunforward(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Reset the domain to a homogeneous medium for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[node,face,elem]=meshabox([40 0 0], [160, 120, 60], 10);
cfg=rbsetmesh(cfg,node,elem,cfg.prop,ones(size(node,1),1));

% [nosp,fcsp]=meshasphere(s0, 5, 3);
% [no,fc]=mergemesh(nobbx, fcbbx, nosp, fcsp);
% 
% [cfg.node, cfg.elem]=s2m(no,fc(:,1:3),1,40,'tetgen',[41 1 1;s0]);
% cfg=rbmeshprep(cfg);

sd=rbsdmap(cfg);

[recon.node,face,recon.elem]=meshabox([40 0 0], [160, 120, 60], 40);

[f2rid, f2rweight]=tsearchn(recon.node,recon.elem,cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run 10 iterations to recover mua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxiter=10;
resid=zeros(1,maxiter);

% initialize reconstruction to homogeneous
recon.prop=cfg.prop(ones(size(recon.node,1),1),:);

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
    recon.prop(:,1)=recon.prop(:,1) + dmu_recon(:);          % update forward mesh mua vector
    cfg.prop=meshinterp(recon.prop,f2rid, f2rweight,recon.elem); % interpolate the update to the forward mesh
    fprintf(1,'iter [%4d]: residual=%e, relres=%e (time=%f s)\n',i, resid(i), resid(i)/resid(1), toc);
end

plotmesh([cfg.node,cfg.mua],cfg.elem,'z=20','facecolor','interp','linestyle','none')
hold on;
plotmesh([cfg.node,cfg.mua],cfg.elem,'x=70','facecolor','interp','linestyle','none')
view(3);