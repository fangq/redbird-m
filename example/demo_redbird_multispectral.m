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

cfg0.prop = containers.Map();
cfg0.prop('690')=[
    0 0 1 1
    0.008 1 0 1.37
    0.016 1 0 1.37
];
cfg0.prop('830')=[
    0 0 1 1
    0.008 1 0 1.37
    0.016 1 0 1.37
];

wavelengths=cfg0.prop.keys;

%z0=1/(cfg0.prop(2,1)+cfg0.prop(2,2)*(1-cfg0.prop(2,3)));
%cfg0.srcpos(:,3)=cfg0.srcpos(:,3)+z0;
%cfg0.detpos(:,3)=cfg0.detpos(:,3)-z0;

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

cfg.mua=containers.Map();

for waveid=wavelengths
   wv=waveid{1};
   prop=cfg.prop(wv);
   cfg.mua(wv)=ones(size(cfg.node,1),1)*prop(cfg.seg(1)+1,1);
end

for i=1:maxiter
    tic
    [dmu, misfit]=rbreconstep(cfg,sd,recon,detphi0,f2rid,f2rweight);
    resid(i)=sum(abs(misfit));         % store the residual
    cfg.mua=cfg.mua + dmu(:);          % update forward mesh mua vector
    fprintf(1,'iter [%4d]: residual=%e, relres=%e (time=%f s)\n',i, resid(i), resid(i)/resid(1), toc);
end

plotmesh([cfg.node,cfg.mua],cfg.elem,'z=20','facecolor','interp','linestyle','none')
hold on;
plotmesh([cfg.node,cfg.mua],cfg.elem,'x=70','facecolor','interp','linestyle','none')
view(3);