%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography,
%      Copyright Qianqina Fang, 2018
%
% In this example, we show the most basic usage of Redbird.
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('rbrun', 'file'))
    addpath(fullfile(pwd, '../matlab'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg xcfg;

[cfg.node, cfg.elem] = meshgrid6(0:2:60, 0:2:60, 0:2:30);
cfg.elem(:, 1:4) = meshreorient(cfg.node(:, 1:3), cfg.elem(:, 1:4));
cfg.face = volface(cfg.elem);

nn = size(cfg.node, 1);
cfg.seg = ones(size(cfg.elem, 1), 1);

cfg.srctype = 'planar';
cfg.srcpos = [9.5 9.5 0];
cfg.srcparam1 = [40 0 0 0];
cfg.srcparam2 = [0 40 0 0];
cfg.srcdir = [0 0 1];

cfg.dettype = 'planar';
cfg.detpos = [10 10 30];
cfg.detparam1 = [40 0 0 0];
cfg.detparam2 = [0 40 0 0];
cfg.detdir = [0 0 -1];

cfg.prop = [0 0 1 1; 0.005 1 0 1.37];
cfg.omega = 0;
cfg.srcweight = 2;

z0 = 1 / (cfg.prop(2, 1) + cfg.prop(2, 2) * (1 - cfg.prop(2, 3)));

cfg = rbmeshprep(cfg);

save config.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build LHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Amat, deldotdel] = rbfemlhs(cfg, cfg.deldotdel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rhs, loc, bary] = rbfemrhs(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Solve for solutions at all freenodes: Afree*sol=rhs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf(1, 'solving for the solution ...\n');
% phi=rbfemsolve(Amat,rhs,'qmr',1e-6,100);
phi = rbfemsolve(Amat, rhs);
phi(phi < 0) = 0;
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Extract detector readings from the solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detval = rbfemgetdet(phi, cfg, rhs); % or detval=rbfemgetdet(phi, cfg, loc, bary);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('mcxlab', 'file'))
    xcfg.nphoton = 1e8;
    xcfg.vol = uint8(ones(60, 60, 30));
    xcfg.srcdir = [0 0 1 0];
    xcfg.gpuid = 1;
    xcfg.autopilot = 1;
    xcfg.prop = cfg.prop;
    xcfg.tstart = 0;
    xcfg.tend = 5e-9;
    xcfg.tstep = 5e-9;
    xcfg.seed = 99999;
    xcfg.issrcfrom0 = 0;

    % a uniform planar source outside the volume
    xcfg.srctype = 'planar';
    xcfg.srcpos = [10 10 0];
    xcfg.srcparam1 = [40 0 0 0];
    xcfg.srcparam2 = [0 40 0 0];

    flux = mcxlab(xcfg);
    fcw = flux.data * xcfg.tstep;
    subplot(211);
    imagesc(rot90(log10(abs(squeeze(cfg.srcweight * fcw(:, 30, :))))));
    axis equal;
    colorbar;
    set(gca, 'xlim', [0 60], 'ylim', [0 30]);
    title('MCX solution');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cl = get(subplot(211), 'clim');
subplot(212);
plotmesh([cfg.node full(log10(phi(:, 1)))], cfg.elem, 'x>30');
view([-1 0 0]);
shading interp;
set(gca, 'clim', cl);
colorbar;
title('Redbird solution');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clines = 0.5:-0.5:-5;
[xi, yi] = meshgrid(0.5:59.5, 0.5:29.5);
[cutpos, cutvalue, facedata] = qmeshcut(cfg.elem, cfg.node, phi(:, 1), 'x=29.5');
vphi = griddata(cutpos(:, 2), cutpos(:, 3), cutvalue, xi + 0.5, yi);

figure;
[c, h] = contour(xi, yi, log10(vphi), clines, 'r-', 'LineWidth', 2);

cwf = squeeze(fcw(30, :, :))';
hold on;
contour(xi, yi, log10(cfg.srcweight * cwf), clines, 'b-', 'LineWidth', 2);
