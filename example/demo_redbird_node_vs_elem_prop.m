clear cfg

[cfg.node,cfg.face, cfg.elem]=meshabox([40 0 0], [160, 120, 60], 10);

nn=size(cfg.node,1);

[xi,yi]=meshgrid(60:20:140,20:20:100);
cfg.srcpos=[xi(:),yi(:),zeros(numel(yi),1)];
cfg.detpos=[xi(:),yi(:),60*ones(numel(yi),1)];
cfg.srcdir=[0 0 1];
cfg.detdir=[0 0 -1];

% label based properties for elements
cfg.prop=[
    0 0 1 1
    0.008 1 0 1.37
];
cfg.seg=ones(size(cfg.elem,1),1);

cfg.omega=2*pi*70e6;

tic
cfg=rbmeshprep(cfg);
fprintf('preparing mesh ... \t%f seconds\n',toc);

rbgetbulk(cfg)

%% RF forward using segmentation based properties

[detphi1,phi1]=rbrunforward(cfg);

%%  node-based properties

% cfg.seg=ones(size(cfg.node,1),1);
% cfg.prop=[
%     0 0 1 1
%     0.008 1 0 1.37
%     0.016 1 0 1.37
% ];
cfg.prop=repmat([0.008 1 0 1.37],size(cfg.node,1),1);

rbgetbulk(cfg)

[detphi2,phi2]=rbrunforward(cfg);

%% comparing the two solutions

rr=detphi1./detphi2;
fprintf(1,'max ratio=%f , min ratio=%f\n',max(abs(rr(:))), min(abs(rr(:))));
fprintf(1,'max angle diff=%f , min angle diff=%f\n',max(angle(rr(:))), min(angle(rr(:))));
