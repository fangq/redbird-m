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

clear cfg xcfg

[cfg.node, cfg.elem]=meshgrid5(0:2:60,0:2:60,0:2:30);
cfg.face=volface(cfg.elem);

nn=size(cfg.node,1);
cfg.seg=ones(size(cfg.elem,1),1);
c0=meshcentroid(cfg.node,cfg.face);

idx=find(c0(:,3)==0 & c0(:,1)>10 & c0(:,1)<50 & c0(:,2)>10 & c0(:,2)<50);

cfg.srcpos=zeros(1,size(cfg.face,1));
cfg.srcpos(idx)=1;

cfg.srcdir=[0 0 1];

[xi,yi]=meshgrid(10:10:50,20:10:40);
cfg.detpos=[xi(:),yi(:),30*ones(numel(yi),1)];

cfg.prop=[0 0 1 1;0.005 1 0 1.37];
cfg.omega=0;

cfg=rbmeshprep(cfg);

save config.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build LHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Amat,deldotdel]=rbfemlhs(cfg);
%[deldotdel2]=rbdeldotdel(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rhs,loc,bary]=rbfemrhs(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Solve for solutions at all freenodes: Afree*sol=rhs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;fprintf(1,'solving for the solution ...\n');
%phi=rbfemsolve(Amat,rhs,'qmr',1e-6,100);
phi=rbfemsolve(Amat,rhs);
phi(phi<0)=0;
toc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Extract detector readings from the solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%detval=rbfemgetdet(phi, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs); 

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

if(exist('mcxlab','file'))
        xcfg.nphoton=1e8;
        xcfg.vol=uint8(ones(60,60,30));
        xcfg.srcdir=[0 0 1 0];
        xcfg.gpuid=2;
        xcfg.autopilot=1;
        xcfg.prop=cfg.prop;
        xcfg.tstart=0;
        xcfg.seed=99999;

        % a uniform planar source outside the volume
        xcfg.srctype='planar';
        xcfg.srcpos=[10 10 0];
        xcfg.srcparam1=[40 0 0 0];
        xcfg.srcparam2=[0 40 0 0];
        xcfg.tend=5e-9;
        xcfg.tstep=5e-9;
        flux=mcxlab(xcfg);
        fcw=flux.data*xcfg.tstep;
        subplot(211);
        imagesc(rot90(log10(abs(squeeze(fcw(:,30,:))))))
        axis equal; colorbar
        set(gca,'xlim',[0 60],'ylim',[0 30])
        title('MCX solution');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if you have the SVN version of iso2mesh, use the next line to plot:
% qmeshcut(cfg.elem(:,1:4),cfg.node(:,1:3),log10(abs(flux.data(:))),'y=30','linestyle','none');

cl=get(subplot(211),'clim');
subplot(212);
plotmesh([cfg.node full(log10(phi(:)))],cfg.elem,'x>30')
view([-1 0 0]);
shading interp;
set(gca,'clim',cl);
colorbar;
title('Redbird solution');