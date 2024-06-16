%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography, 
%      Copyright Qianqina Fang, 2018
%
% Continuous-Wave (CW) reconstruction of absorption (mua) target
% (streamlined version by calling rbrunrecon)
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('rbrun','file'))
    addpath(fullfile(pwd, '../matlab'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg cfg0 recon

s0=[100, 50, 30];

[nobbx,fcbbx]=meshabox([40 0 0], [160, 120, 40], 10);
[nosp,fcsp]=meshasphere(s0, 5, 1);
[no,fc]=mergemesh(nobbx, fcbbx, nosp, fcsp);

[cfg0.node, cfg0.elem]=s2m(no,fc(:,1:3),1,40,'tetgen1.5',[41 1 1;s0]);
nn=size(cfg0.node,1);
cfg0.seg=cfg0.elem(:,5);
cfg0.srcdir=[0 0 1];

[xi,yi]=meshgrid(50:10:90,30:10:70);
cfg0.srcpos=[xi(:),yi(:),zeros(numel(yi),1)];
cfg0.detpos=[xi(:),yi(:),40*ones(numel(yi),1)];
cfg0.detdir=[0 0 -1];

cfg0.prop=[
    0 0 1 1
    0.008 1 0 1.37
    0.016 1 0 1.37
];

cfg0.omega=0;%67.5e6*2*pi;

cfg=cfg0;

cfg0=rbmeshprep(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run forward for the heterogeneous domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detphi0=rbrun(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Reset the domain to a homogeneous medium for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create forward mesh for reconstruction - if it has a different density
% compared to the forward mesh used for generating data, it may introduce
% numerical error - so, set density to 10 gives the best result, need
% to debug this further

[node,face,elem]=meshabox([40 0 0], [160, 120, 40], 10);
cfg=rbsetmesh(cfg,cfg0.node,cfg0.elem,cfg.prop,ones(size(cfg0.node,1),1));

sd=rbsdmap(cfg);

% create coarse reconstruction mesh
[recon.node,face,recon.elem]=meshabox([40 0 0], [160, 120, 40], 20);
[recon.mapid, recon.mapweight]=tsearchn(recon.node,recon.elem,cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Streamlined reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize reconstruction to homogeneous (label=1)
recon.prop=cfg.prop(ones(size(recon.node,1),1)+1,:);
cfg.prop=cfg.prop(ones(size(cfg.node,1),1)+1,:);
cfg=rmfield(cfg,'seg');

% run stream-lined image reconstruction
[newrecon,resid,newcfg,~,~,phi]=rbrun(cfg,recon,detphi0,sd,'mex',0,'maxiter',10,'lambda',1e-4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
plotmesh([newcfg.node,newcfg.prop(:,1)],newcfg.elem,'z=30','facecolor','interp','linestyle','none')
colorbar
hold on;
plotmesh([newcfg.node,newcfg.prop(:,1)],newcfg.elem,'x=70','facecolor','interp','linestyle','none')
view(3);
