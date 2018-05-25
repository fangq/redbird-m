%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography, 
%      Copyright Qianqina Fang, 2018
%
% In this example, we show the most basic usage of Redbird.
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../matlab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg

% [cfg.node,cfg.elem]=readrbmesh('lightphantom_newf');
% cfg.node=cfg.node*1000;

[cfg.node,cfg.face, cfg.elem]=meshabox([40 0 0], [160, 120, 60], 10);

%[cfg.node, cfg.face, cfg.elem]=meshabox([0 0 0],[60 60 30],3);
nn=size(cfg.node,1);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcdir=[0 0 1];

[xi,yi]=meshgrid(60:20:140,20:20:100);
cfg.srcpos=[xi(:),yi(:),zeros(numel(yi),1)];
cfg.detpos=[xi(:),yi(:),60*ones(numel(yi),1)];

cfg.prop=[
    0 0 1 1
    0.008 1 0 1.37
    0.016 1 0 1.37
];

z0=1/(cfg.prop(2,1)+cfg.prop(2,2)*(1-cfg.prop(2,3)));

cfg.srcpos(:,3)=cfg.srcpos(:,3)+z0;
cfg.detpos(:,3)=cfg.detpos(:,3)-z0;

cfg.omega=2*pi*70e6;
cfg.omega=0;

cfg=rbmeshprep(cfg);

save config.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build LHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% [Amat,deldotdel]=rbfemlhs(cfg); % use mex function rbfemmatrix, 5x faster
% toc

tic
deldotdel=rbdeldotdel(cfg);
Amat=rbfemlhs(cfg,deldotdel); % use native matlab code, 1 sec for 50k nodes
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rhs,loc,bary]=rbfemrhs(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Solve for solutions at all freenodes: Afree*sol=rhs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;fprintf(1,'solving for the solution ...\n');
%phi=rbfemsolve(Amat,rhs,'symmlq',1e-20,100);
phi=rbfemsolve(Amat,rhs);
toc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Extract detector readings from the solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detval=rbfemgetdet(phi, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build Jacobians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nvol=nodevolume(cfg.node,cfg.elem);
% sd=rbsdmap(cfg);
% Jmua=rbjacmua(sd, phi, nvol);
% Jd=rbjacdcoef(sd, phi, deldotdel, cfg.elem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sid=13;

cfg.srcpos(:,3)=cfg.srcpos(:,3)-z0;

srcloc=cfg.srcpos(sid,1:3);
detloc=cfg.node;
phicw=cwdiffusion(cfg.prop(2,1), cfg.prop(2,2)*(1-cfg.prop(2,3)), 0.431, srcloc, detloc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(221);plotmesh([cfg.node,log10(abs(phi(1:size(cfg.node,1),sid)))],cfg.elem,'y=30','facecolor','interp','linestyle','none')
cl=get(gca,'clim');
set(gca, 'xlim', [60, 140]);
set(gca, 'zlim', [0 60]);
view([0 1 0]);
colorbar;

subplot(222);plotmesh([cfg.node,log10(abs(phicw(1:size(cfg.node,1),1)))],cfg.elem,'y=30','facecolor','interp','linestyle','none')
view([0 1 0]);
set(gca, 'xlim', [60, 140]);
set(gca, 'zlim', [0 60]);
set(gca, 'clim', cl);
colorbar;

dd=log10(abs(phi(1:size(cfg.node,1),sid))) - log10(abs(phicw(1:size(cfg.node,1),1)));
subplot(223);plotmesh([cfg.node,dd],cfg.elem,'y=30','facecolor','interp','linestyle','none')
view([0 1 0]);
set(gca, 'xlim', [60, 140]);
set(gca, 'zlim', [0 60]);
colorbar;

subplot(224);plotmesh([cfg.node,dd],cfg.elem,'y=30','facecolor','interp','linestyle','none')
hist(dd(:),100);