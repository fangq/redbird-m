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

clear cfg cfg0 recon;

s0 = [70, 50, 20];

[nobbx, fcbbx] = meshabox([40 0 0], [160, 120, 60], 30);
[nosp, fcsp] = meshasphere(s0, 5, 1);
[no, fc] = mergemesh(nobbx, fcbbx, nosp, fcsp);

[cfg0.node, cfg0.elem] = s2m(no, fc(:, 1:3), 1, 40, 'tetgen', [41 1 1; s0]);

cfg0.seg = cfg0.elem(:, 5);
cfg0.srcdir = [0 0 1];

[xi, yi] = meshgrid(60:20:140, 20:20:100);
cfg0.srcpos = [xi(:), yi(:), zeros(numel(yi), 1)];
cfg0.detpos = [xi(:), yi(:), 60 * ones(numel(yi), 1)];
cfg0.detdir = [0 0 -1];

cfg0.prop = [
             0 0 1 1
             0.008 1 0 1.37
             0.016 1 0 1.37
            ];

cfg0.omega = 2 * pi * 70e6;

cfg = cfg0;

cfg0 = rbmeshprep(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run forward for the heterogeneous domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detphi0 = rbrunforward(cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Reset the domain to a homogeneous medium for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create forward mesh for reconstruction - if it has a different density
% compared to the forward mesh used for generating data, it may introduce
% numerical error - so, set density to 10 gives the best result, need
% to debug this further

[node, face, elem] = meshabox([40 0 0], [160, 120, 60], 30);
cfg = rbsetmesh(cfg, node, elem, cfg.prop, ones(size(node, 1), 1));

sd = rbsdmap(cfg);

% create coarse reconstruction mesh
[recon.node, face, recon.elem] = meshabox([40 0 0], [160, 120, 60], 40);
[recon.mapid, recon.mapweight] = tsearchn(recon.node, recon.elem, cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Bulk fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize reconstruction to homogeneous (label=1)
recon.bulk = struct('mua', 0.003, 'musp', 0.8);

% calling rbrunrecon is equivalent to calling the below for-loop
newrecon = rbrun(cfg, recon, detphi0, sd, 'bulk');

newrecon.prop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Image reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize reconstruction to homogeneous (label=1)
recon.bulk = struct('mua', newrecon.prop(2, 1), 'musp', newrecon.prop(2, 2));
% recon.bulk=struct('mua',0.008, 'musp', 1);

recon.prop = [];

% calling rbrunrecon is equivalent to calling the below for-loop
[newrecon, resid, newcfg] = rbrun(cfg, recon, detphi0, sd, 'image');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotmesh([newcfg.node, newcfg.prop(:, 1)], newcfg.elem, 'z=20', 'facecolor', 'interp', 'linestyle', 'none');
hold on;
plotmesh([newcfg.node, newcfg.prop(:, 1)], newcfg.elem, 'x=70', 'facecolor', 'interp', 'linestyle', 'none');
view(3);
