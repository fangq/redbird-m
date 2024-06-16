if(~exist('rbrun','file'))
    addpath(fullfile(pwd, '../matlab'));
end

%% Generate Source/Detector Patterns

srcpattern = omciBarPattern(16, 16, 1);
detpattern=srcpattern;

%%
clear cfg cfg0 recon

% Centers of two spherical inclusions
s0=[40 40 20];
s2=[90 20 20];

% Bounding box w/ 2 inclusions
[nobbx,fcbbx]=meshabox([0 -10 0], [120 70 40], 10);
[nosp,fcsp]=meshasphere(s0, 7.5, 1);
[nosp2,fcsp2] = meshasphere(s2,7.5,1);
[no,fc]=mergemesh(nobbx, fcbbx, nosp, fcsp);
[no,fc] = mergemesh(no,fc,nosp2,fcsp2);

[cfg0.node, cfg0.elem]=s2m(no,fc(:,1:3),1,10,'tetgen',[11 1 1;s0;s2]);

cfg0.srctype = 'pattern';   % Source setup
cfg0.srcpos = [10 10 0];
cfg0.srcparam1 = [100 0 0 0];
cfg0.srcparam2 = [0 40 0 0];
cfg0.srcdir = [0 0 1];
cfg0.srcpattern = srcpattern; % # src of patterns based upon  size(srcpattern,3)

cfg0.dettype = 'pattern';   % Detector setup
cfg0.detpos = [10 10 40];
cfg0.detparam1 = [100 0 0 0];
cfg0.detparam2 = [0 40 0 0];
cfg0.detdir = [0 0 -1];
cfg0.detpattern = srcpattern;   % Same patterns used as sources

cfg0.prop=[
    0 0 1 1
    0.008 1 0 1.37
    0.016 1 0 1.37  % Currently using same OP for both inclusions
    0.008 2 0 1.37
];

cfg0.omega=2*pi*65e6;
% cfg0.omega=0;   % CW (no frequency modulation)

cfg=cfg0;

cfg0=rbmeshprep(cfg0);

%% run forward simulation to get simulated data

[detphi0,phi0] =rbrun(cfg0);

%% create reconstruction dual-mesh

[node,face,elem] = meshabox([0 -10 0], [120 70 40], 10);
cfg = rbsetmesh(cfg,cfg.node,cfg.elem,[0 0 1 1;0.010 1 0 1.37],ones(size(cfg.node,1),1));
sd = rbsdmap(cfg);    % Acquire sd map for calculating Jmua

% Create a coarser reconstruction mesh for faster inversion
[recon.node,face,recon.elem] = meshabox([0 -10 0], [120 70 40], 30);
[recon.mapid, recon.mapweight] = tsearchn(recon.node,recon.elem,cfg.node);
recon.bulk=struct('mua',0.008,'musp',1);

%% run image reconstruction

[newrecon,resid] = rbrun(cfg,recon,detphi0,'lambda',1e-6,'maxiter',5);


%% plot reconstructed images
figure
plotmesh([newrecon.node,newrecon.prop(:,1)],newrecon.elem,'z=20','facecolor','interp','linestyle','none')
hold on;
plotmesh([newrecon.node,newrecon.prop(:,1)],newrecon.elem,'x=90','facecolor','interp','linestyle','none')
plotmesh([newrecon.node,newrecon.prop(:,1)],newrecon.elem,'x=40','facecolor','interp','linestyle','none')
title('Reconstructed Absorption (\mu_{a})');
set(gca, 'FontSize', 16, 'FontWeight', 'bold');
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
hold off;
view(3);

figure
plotmesh([newrecon.node,newrecon.prop(:,2)],newrecon.elem,'z=20','facecolor','interp','linestyle','none')
hold on;
plotmesh([newrecon.node,newrecon.prop(:,2)],newrecon.elem,'x=90','facecolor','interp','linestyle','none')
plotmesh([newrecon.node,newrecon.prop(:,2)],newrecon.elem,'x=40','facecolor','interp','linestyle','none')
title('Reconstructed Scattering (\mu_{s}'')');
set(gca, 'FontSize', 16, 'FontWeight', 'bold');
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
hold off; 
view(3);

%%

[xi,yi] = meshgrid(0.5:119.5,-9.5:69.5);
[cutpos,cutvalue,facedata] = qmeshcut(newrecon.elem,newrecon.node,newrecon.prop(:,1),'z=20');
vphi = griddata(cutpos(:,1),cutpos(:,2),cutvalue,xi+0.5,yi);

figure,imagesc(vphi)