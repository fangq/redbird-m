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

srcpattern = diag(ones(1, 8));
srcpattern = permute(repmat(srcpattern, [1, 1, size(srcpattern, 1)]), [2 3 1]);
srcpattern = cat(3, srcpattern, permute(srcpattern, [2 1 3]));
detpattern = srcpattern;

s0 = [80, 40, 20]; % center of the inclusion (in mm)
rs = 5;            % radius of the sphere (in mm)

[nobbx, fcbbx] = meshabox([40 0 0], [140 100 50], 25);
[nosp, fcsp] = meshasphere(s0, rs, 1);
[no, fc] = mergemesh(nobbx, fcbbx, nosp, fcsp);

[cfg0.node, cfg0.elem] = s2m(no, fc(:, 1:3), 1, 25, 'tetgen', [41 1 1; s0]);

% [cfg.node, cfg.face, cfg.elem]=meshabox([0 0 0],[60 60 30],3);
nn = size(cfg0.node, 1);
cfg0.seg = cfg0.elem(:, 5);
cfg0.srcdir = [0 0 1];

[xi, yi] = meshgrid(60:10:120, 20:10:80);
cfg0.srcpos = [xi(:), yi(:), zeros(numel(yi), 1)];
cfg0.detpos = [xi(:), yi(:), 50 * ones(numel(yi), 1)];
cfg0.detdir = [0 0 -1];

cfg0.srcpos = [cfg0.srcpos; 60 20 0];
cfg0.srcid = length(xi(:)) + 1;
cfg0.srcparam1 = [60 0 0 0];
cfg0.srcparam2 = [0 60 0 0];
cfg0.srctype = 'pattern';
cfg0.srcpattern = srcpattern;
cfg0.srcweight = ones(1, 16);

cfg0.detpos = [cfg0.detpos; 60 20 50];
cfg0.detid = length(xi(:)) + 1;
cfg0.detparam1 = [60 0 0 0];
cfg0.detparam2 = [0 60 0 0];
cfg0.dettype = 'pattern';
cfg0.detpattern = srcpattern;
cfg0.detweight = ones(1, 16);

cfg0.param = struct;
cfg0.param.hbo = [10 20];
cfg0.param.hbr = [4  8];
cfg0.param.scatamp = [1.6e-6 2.4e-6];
cfg0.param.scatpow = [0.9177 0.9177];

cfg0.prop = containers.Map();  % if both prop and param are defined, param will ovewrite prop
cfg0.prop('690') = [0 0 1 1; 0.006 0.8 0 1.37; 0.012 1 0 1.37];
cfg0.prop('830') = [0 0 1 1; 0.005 0.5 0 1.37; 0.010 0.6 0 1.37];

cfg0.wavesrc = containers.Map({'690', '830'}, {[1:length(xi(:)) + 1], [1:length(xi(:)) + 1]});
cfg0.rfcw.src = containers.Map({'RF', 'CW'}, {[1:length(xi(:))], [length(xi(:)) + 1]});

cfg0.wavedet = containers.Map({'690', '830'}, {[1:length(xi(:)) + 1], [1:length(xi(:)) + 1]});
cfg0.rfcw.det = containers.Map({'RF', 'CW'}, {[1:length(xi(:))], [length(xi(:)) + 1]});

wavelengths = cfg0.prop.keys;

cfg0.omega = containers.Map({'690', '830'}, {67.5e6 * 2 * pi, 2 * pi * 75e6});

cfg = cfg0;

%

[cfg0, sd0] = rbmeshprep(cfg0);
rfcw = [1 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   run forward for all wavelengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[detphi0, phi0] = rbrunforward(cfg0, 'sd', sd0, 'rfcw', rfcw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   run reconstruction using the forward data, setup dual-mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[node, face, elem] = meshabox([40 0 0], [140, 100, 50], 25);
clear face;

cfg = rbsetmesh(cfg, node, elem, cfg.prop, ones(size(node, 1), 1));

[recon.node, face, recon.elem] = meshabox([40 0 0], [140, 100, 50], 50);
clear face;
[recon.mapid, recon.mapweight] = tsearchn(recon.node, recon.elem, cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   run bulk fitting first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sd = rbsdmap(cfg);
recon.bulk = struct('hbo', 8, 'hbr', 2, 'scatamp', 1.8e-6, 'scatpow', 0.9); % Required: this gives initial guesses
recon.param = struct('hbo', 8, 'hbr', 2, 'scatamp', 1.8e-6, 'scatpow', 0.9); % Required: this defines chromophores
recon.prop = containers.Map({'690', '830'}, {[], []}); % Required: for wavelengths

%%

[newrecon, resid] = rbrun(cfg, recon, detphi0, sd, 'mode', 'bulk', 'rfcw', rfcw, 'lambda', 1e-3, 'maxiter', 5);

recon.bulk = newrecon.param;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   take the fitted bulk and set it for full image recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[newrecon, resid, newcfg] = rbrun(cfg, recon, detphi0, sd, 'mode', 'image', 'rfcw', rfcw, 'lambda', 1e-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plotmesh([newrecon.node, newrecon.param.hbo(:) + newrecon.param.hbr(:)], newrecon.elem, 'z=20', 'facecolor', 'interp', 'linestyle', 'none');
hold on;
plotmesh([newrecon.node, newrecon.param.hbo(:) + newrecon.param.hbr(:)], newrecon.elem, 'x=90', 'facecolor', 'interp', 'linestyle', 'none');
view(3);
