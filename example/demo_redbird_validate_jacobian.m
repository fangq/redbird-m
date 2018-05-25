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

[cfg.node,cfg.face, cfg.elem]=meshabox([40 0 0], [160, 120, 60], 20);

%[cfg.node, cfg.face, cfg.elem]=meshabox([0 0 0],[60 60 30],3);
nn=size(cfg.node,1);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcdir=[0 0 1];

[xi,yi]=meshgrid(60:20:140,20:20:100);
cfg.srcpos=[xi(:),yi(:),zeros(numel(yi),1)];
cfg.detpos=[xi(:),yi(:),60*ones(numel(yi),1)];
cfg.detpos(1:5,:)=[];

cfg.prop=[
    0 0 1 1
    0.008 1 0 1.37
    0.00801 1 0 1.37
];

sid=13;
zdepth=8;
srad=2;

c0=meshcentroid(cfg.node,cfg.elem);

perturbeid=find(abs(c0(:,1)-cfg.srcpos(sid,1))<srad & abs(c0(:,2)-cfg.srcpos(sid,1))<srad & abs(c0(:,3)-zdepth)<srad);
perturbeid=perturbeid(1);

z0=1/(cfg.prop(2,1)+cfg.prop(2,2)*(1-cfg.prop(2,3)));

cfg.srcpos(:,3)=cfg.srcpos(:,3)+z0;
cfg.detpos(:,3)=cfg.detpos(:,3)-z0;

cfg.omega=2*pi*70e6;
cfg.omega=0;

cfg=rbmeshprep(cfg);

%save config.mat

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
%%   perturb mua in an element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.elemprop(perturbeid)=2;

tic
Amat2=rbfemlhs(cfg,deldotdel); % use native matlab code, 1 sec for 50k nodes
toc

phi2=rbfemsolve(Amat2,rhs);
detval2=rbfemgetdet(phi2, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   perturb D in an element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.elemprop(perturbeid)=2;
cfg.prop(3,:)=cfg.prop(2,:);
cfg.prop(3,2)=cfg.prop(3,2)*1.01;

tic
Amat3=rbfemlhs(cfg,deldotdel); % use native matlab code, 1 sec for 50k nodes
toc

phi3=rbfemsolve(Amat3,rhs);
detval3=rbfemgetdet(phi3, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build mua Jacobians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isnodal=1;  % isnodal=1 builds nodal Jacobain; isnodal=0 builds elem-based Jacobians

nvol=nodevolume(cfg.node,cfg.elem, cfg.evol);
sd=rbsdmap(cfg);

if(isnodal)
    Jmua=rbjacmua(sd, phi, nvol); % build nodal-based Jacobian for mua
else
    Jmua=rbjacmua(sd, phi, cfg.evol, cfg.elem); % build elem-based J_mua, large & slow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Compare mua Jacobian with direct measurement change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dmua=(cfg.prop(3,1)-cfg.prop(2,1));
dphi1=(detval2-detval);          % change of measurement from separate forward
if(isnodal)
    totalvol=sum(nvol(cfg.elem(perturbeid,:)));
    dphi2=sum(Jmua(:,cfg.elem(perturbeid,:)),2)*dmua*(cfg.evol(perturbeid)/totalvol);
else
    dphi2=Jmua(:,perturbeid)*dmua;  % change of measurement predicted from Jacobians
end
dphi2=reshape(dphi2,size(cfg.detpos,1),size(cfg.srcpos,1));

figure;
plot(1:length(dphi1),dphi1,'r-o',1:length(dphi2),dphi2,'b-+');
legend('direct measurement change','predicted from J_{\mua}');

figure;
dd=dphi1./dphi2;
fprintf(1,'isnodal=\t%d\tsum(Jmua)=\t%f\tmedian(dphi1/dphi2)=\t%f\n',isnodal,sum(Jmua(:)), median(dd(:)));
imagesc(dd)
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build D Jacobians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[Jd, JD]=rbjacdcoef(sd, phi, deldotdel, cfg.elem); % build nodal-based Jacobian for mua
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Compare D Jacobian with direct measurement change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isnodal=0;

dD=(1/3/cfg.prop(3,2)-1/3/cfg.prop(2,2));
dphi3=(detval3-detval);          % change of measurement from separate forward
if(isnodal)
    totalvol=sum(nvol(cfg.elem(perturbeid,:)));
    dphi4=sum(Jd(:,cfg.elem(perturbeid,:)),2)*dD*(cfg.evol(perturbeid)/totalvol);
else
    dphi4=JD(:,perturbeid)*dD;  % change of measurement predicted from Jacobians
end
dphi4=reshape(dphi4,size(cfg.detpos,1),size(cfg.srcpos,1));

figure;
plot(1:length(dphi3),dphi3,'r-o',1:length(dphi4),dphi4,'b-+');
legend('direct measurement change','predicted from J_{D}');

figure;
dd2=dphi3./dphi4;
fprintf(1,'isnodal=\t%d\tsum(Jd)=\t%f\tmedian(dphi3/dphi4)=\t%f\n',isnodal,sum(Jd(:)), median(dd2(:)));
imagesc(dd2)
colorbar;