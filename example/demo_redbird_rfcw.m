%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redbird - A Diffusion Solver for Diffuse Optical Tomography, 
%      Copyright Qianqina Fang, 2018
%
% In this example, we show the most basic usage of Redbird.
%
% This file is part of Redbird URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if(~exist('rbrun','file'))
%     addpath(fullfile(pwd, '../matlab'));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepare simulation input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cfg cfg0 recon

srcpattern = diag(ones(1,16))+diag(ones(1,15),-1);
srcpattern(1,end)=1;
srcpattern=permute(repmat(srcpattern,[1,1,16]),[2 3 1]);
srcpattern=cat(3,srcpattern,permute(srcpattern,[2 1 3]));
detpattern=srcpattern;

s0=[90, 50, 20]; % center of the inclusion (in mm)
rs=5;            % radius of the sphere (in mm)

[nobbx,fcbbx]=meshabox([40 0 0], [160, 120, 60], 10);
[nosp,fcsp]=meshasphere(s0, rs, 1);
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

cfg0.srcpos = [cfg0.srcpos; 60 20 0];
cfg0.srcid = 26;
cfg0.srcparam1 = [80 0 0 0];
cfg0.srcparam2 = [0 80 0 0];
cfg0.srctype = 'pattern';
cfg0.srcpattern = srcpattern;
cfg0.srcweight = rand(1,32);

cfg0.detpos = [cfg0.detpos; 60 20 60];
cfg0.detid = 26;
cfg0.detparam1 = [80 0 0 0];
cfg0.detparam2 = [0 80 0 0];
cfg0.dettype = 'pattern';
cfg0.detpattern = srcpattern;
cfg0.detweight = rand(1,32);
% cfg0.widesrcid = containers.Map();
% cfg0.widesrcid('690') = struct('srcid',26,'srcpattern',srcpattern,'srcweight',srcweight);
% cfg0.widesrcid('830') = struct('srcid',26,'srcpattern',srcpattern,'srcweight',srcweight);

cfg0.param=struct;
cfg0.param.hbo=[15 45];
cfg0.param.hbr=[4  12];
cfg0.param.scatamp = [2.5 4];
cfg0.param.scatpow = [3 3];

cfg0.prop = containers.Map();  % if both prop and param are defined, param will ovewrite prop
cfg0.prop('690')=[0 0 1 1; 0   1 0 1.37; 0 1 0 1.37];
cfg0.prop('830')=[0 0 1 1; 0 0.8 0 1.37; 0 0.8 0 1.37];

cfg0.wavesrc = containers.Map({'690','830'},{[1:13,26],[14:25,26]});
cfg0.rfcw.src = containers.Map({'RF','CW'},{[1:25],[26]});

cfg0.wavedet = containers.Map({'690','830'},{[1:13,26],[14:25,26]});
cfg0.rfcw.det = containers.Map({'RF','CW'},{[1:25],[26]});

wavelengths=cfg0.prop.keys;

cfg0.omega=containers.Map({'690','830'},{67.5e6*2*pi,2*pi*75e6});
% cfg0.omega = 2*pi*75e6;
% cfg0.omega=0;

cfg=cfg0;

%

[cfg0,sd0]=rbmeshprep(cfg0);
rfcw = [1 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   run forward for all wavelengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[detphi0,phi0] = rbrunforward(cfg0,'sd',sd0,'rfcw',rfcw);
% detphi0=rbrunrecon(0,cfg0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   run reconstruction using the forward data, setup dual-mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[node,face,elem]=meshabox([40 0 0], [160, 120, 60], 10);
clear face

cfg=rbsetmesh(cfg,node,elem,cfg.prop,ones(size(node,1),1));

[recon.node,face,recon.elem]=meshabox([40 0 0], [160, 120, 60], 40);
clear face
[recon.mapid, recon.mapweight]=tsearchn(recon.node,recon.elem,cfg.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   run bulk fitting first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sd = rbsdmap(cfg);
% sd=rbsdmap(cfg,'wavesrc',containers.Map({'690','830'},{[1:2:25],[1:2:25]}));
recon.bulk=struct('hbo',8,'hbr',2,'scatamp',2.7286,'scatpow',2.2781); % Required: this gives initial guesses
recon.param=struct('hbo',8,'hbr',2,'scatamp',2.7286,'scatpow',2.2781); % Required: this defines chromophores
recon.prop=containers.Map({'690','830'},{[],[]}); % Required: for wavelengths

%%

[newrecon,resid]=rbrun(cfg,recon,detphi0,sd,'mode','bulk','rfcw',rfcw,'lambda',1e-3,'maxiter',15);

recon.bulk = newrecon.param;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   take the fitted bulk and set it for full image recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[newrecon,resid,newcfg]=rbrun(cfg,recon,detphi0,sd,'mode','image','reform','logphase');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
plotmesh([newrecon.node,newrecon.param.hbo(:)+newrecon.param.hbr(:)],newrecon.elem,'z=20','facecolor','interp','linestyle','none')
hold on;
plotmesh([newrecon.node,newrecon.param.hbo(:)+newrecon.param.hbr(:)],newrecon.elem,'x=90','facecolor','interp','linestyle','none')
view(3);