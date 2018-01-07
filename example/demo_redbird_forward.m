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
%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg

[cfg.node, cfg.face, cfg.elem]=meshabox([0 0 0],[60 60 30],6);
nn=size(cfg.node,1);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30 30 1; 10,20,1];
cfg.srcdir=[0 0 1];
cfg.prop=[0 0 1 1;0.005 1 0 1.37];
cfg.omega=2*pi*70e6;

cfg=rbmeshprep(cfg);

save config.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Build LHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Amat,deldotdel]=rbfemlhs(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Build RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhs=rbfemrhs(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Solve for solutions at all freenodes: Afree*sol=rhs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;fprintf(1,'solving for the solution ...\n');
%phi=rbfemsolve(Amat,rhs,'qmr',1e-6,100);
phi=rbfemsolve(Amat,rhs);
toc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if you have the SVN version of iso2mesh, use the next line to plot:
% qmeshcut(cfg.elem(:,1:4),cfg.node(:,1:3),log10(abs(flux.data(:))),'y=30','linestyle','none');

figure;plotmesh([cfg.node(:,1:3),log10(abs(phi(1:size(cfg.node,1),1)))],cfg.elem,'y=30','facecolor','interp','linestyle','none')
view([0 1 0]);
colorbar;
