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

clear cfg;

[cfg.node, cfg.face, cfg.elem] = meshabox([0 0 0], [60 60 30], 1);
nn = size(cfg.node, 1);
cfg.seg = ones(size(cfg.elem, 1), 1);
cfg.srcpos = [30 30 0];
cfg.srcdir = [0 0 1];

cfg.prop = [0 0 1 1; 0.005 1 0 1.37];
cfg.omega = 0;

cfg = rbmeshprep(cfg);

% save config.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build LHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[detphi, phi] = rbrun(cfg);
fprintf('forward solution ... \t%f seconds\n', toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('cwdiffusion', 'file'))
    srcpos = [30 30 0];
    phicw = cwdiffusion(cfg.prop(2, 1), cfg.prop(2, 2) * (1 - cfg.prop(2, 3)), cfg.reff, srcpos, cfg.node);
else
    warning('please download MCX from http://mcx.space to use cwdiffuse.m in mcx/utils');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clines = 0:-0.5:-5;
[xi, yi] = meshgrid(0.5:59.5, 0.5:29.5);
[cutpos, cutvalue] = qmeshcut(cfg.elem, cfg.node, phi(:, 1), 'x=29.5');
vphi = griddata(cutpos(:, 2), cutpos(:, 3), cutvalue, xi + 0.5, yi);

figure;
[c, h] = contour(xi, yi, log10(vphi), clines, 'r-', 'LineWidth', 2);

if (exist('cwdiffusion', 'file'))
    [cutpos, cutvalue] = qmeshcut(cfg.elem, cfg.node, phicw(:), 'x=29.5');
    vphidiffu = griddata(cutpos(:, 2), cutpos(:, 3), cutvalue, xi + 0.5, yi);

    hold on;
    contour(xi, yi, log10(vphidiffu), clines, 'b-', 'LineWidth', 2);

    legend('redbird solution (slab)', 'diffusion (semi-infinite)');
end
