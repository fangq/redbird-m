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
cfg0.detdir=[0 0 -1];

cfg0.param=struct;
cfg0.param.hbo=[15 30];
cfg0.param.hbr=[4  8];

cfg0.prop = containers.Map();  % if both prop and param are defined, param will ovewrite prop
cfg0.prop('690')=[0 0 1 1; 0   1 0 1.37; 0 1 0 1.37];
cfg0.prop('830')=[0 0 1 1; 0 0.8 0 1.37; 0 0.8 0 1.37];

wavelengths=cfg0.prop.keys;

% z0=1/(cfg0.prop(2,1)+cfg0.prop(2,2)*(1-cfg0.prop(2,3)));
% cfg0.srcpos(:,3)=cfg0.srcpos(:,3)+z0;
% cfg0.detpos(:,3)=cfg0.detpos(:,3)-z0;

cfg0.omega=2*pi*70e6;
cfg0.omega=0;

cfg=cfg0;

cfg0=rbmeshprep(cfg0);

%% run forward for all wavelengths
detphi0=rbrunrecon(0,cfg0);

%% run reconstruction using the forward data
[node,face,elem]=meshabox([40 0 0], [160, 120, 60], 10);
clear face

cfg=rbsetmesh(cfg,node,elem,cfg.prop,ones(size(node,1),1));

[recon.node,face,recon.elem]=meshabox([40 0 0], [160, 120, 60], 40);
clear face

recon.param=struct;
recon.param.hbo=8; % initial guess
recon.param.hbr=2;
recon.seg=ones(size(recon.node,1),1);

recon.lambda=0.01;

[recon.mapid, recon.mapweight]=tsearchn(recon.node,recon.elem,cfg.node);

%% run bulk fitting

[newcfg,newrecon]=rbrunrecon(10,cfg,recon,detphi0,rbsdmap(cfg),'lambda',0.0002,'tol',0.001,'report',1);

newrecon.param

%% run image reconstruction

recon.param=struct;
recon.param.hbo=newrecon.param.hbo*ones(size(recon.node,1),1);
recon.param.hbr=newrecon.param.hbr*ones(size(recon.node,1),1);
if(isfield(recon,'seg'))
    recon=rmfield(recon,'seg');
end

cfg.param=struct;
cfg.param.hbo=newrecon.param.hbo*ones(size(cfg.node,1),1);
cfg.param.hbr=newrecon.param.hbr*ones(size(cfg.node,1),1);

[newcfg,newrecon]=rbrunrecon(10,cfg,recon,detphi0,rbsdmap(cfg),'lambda',0.05,'tol',0.01,'report',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotmesh([newrecon.node,newrecon.param.hbo(:)],newrecon.elem,'z=20','facecolor','interp','linestyle','none')
hold on;
plotmesh([newrecon.node,newrecon.param.hbo(:)],newrecon.elem,'x=70','facecolor','interp','linestyle','none')
view(3);