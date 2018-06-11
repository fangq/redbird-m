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

clear cfg

%[cfg.node,cfg.face, cfg.elem]=meshabox([40 0 0], [160, 120, 60], 40, 100);
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
    0.04 1 0 1.37
];

sid=13;
zdepth=30;
srad=2;

c0=meshcentroid(cfg.node,cfg.elem);

perturbeid=find(abs(c0(:,1)-cfg.srcpos(sid,1))<srad & abs(c0(:,2)-cfg.srcpos(sid,1))<srad & abs(c0(:,3)-zdepth)<srad);
perturbeid=perturbeid(1);

z0=1/(cfg.prop(2,1)+cfg.prop(2,2)*(1-cfg.prop(2,3)));

cfg.srcpos(:,3)=cfg.srcpos(:,3)+z0;
cfg.detpos(:,3)=cfg.detpos(:,3)-z0;

%cfg.omega=2*pi*70e6;
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
fprintf('creating deldotdel ... \t%f seconds\n',toc);
Amat=rbfemlhs(cfg,deldotdel); % use native matlab code, 1 sec for 50k nodes
fprintf('creating LHS ... \t%f seconds\n',toc);

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
fprintf('solving forward ... \t%f seconds\n',toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Extract detector readings from the solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detval=rbfemgetdet(phi, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   perturb mua in an element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.elemprop(perturbeid)=2;
cfg.prop(3,:)=cfg.prop(2,:);
cfg.prop(3,1)=cfg.prop(3,1)*1.01;

tic
Amat2=rbfemlhs(cfg,deldotdel); % use native matlab code, 1 sec for 50k nodes
fprintf('solving forward with delta mua ... \t%f seconds\n',toc);

phi2=rbfemsolve(Amat2,rhs);
detval_dmua=rbfemgetdet(phi2, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs); 

dmua=(cfg.prop(3,1)-cfg.prop(2,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   perturb D in an element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.elemprop(perturbeid)=2;
cfg.prop(3,:)=cfg.prop(2,:);
cfg.prop(3,2)=cfg.prop(3,2)*1.01;

dD=(1/3/cfg.prop(3,2)-1/3/cfg.prop(2,2));

tic
Amat3=rbfemlhs(cfg,deldotdel); % use native matlab code, 1 sec for 50k nodes
fprintf('solving forward with delta D ... \t%f seconds\n',toc);

phi3=rbfemsolve(Amat3,rhs);
detval_dd=rbfemgetdet(phi3, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   perturb mua at a node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.elemprop=ones(size(cfg.node,1),1);
cfg.elemprop(cfg.elem(perturbeid,1))=2;
cfg.prop(3,:)=cfg.prop(2,:);
cfg.prop(3,1)=cfg.prop(3,1)*1.01;

tic
Amat2=rbfemlhs(cfg);
fprintf('solving forward with nodal delta mua ... \t%f seconds\n',toc);

phi2_node=rbfemsolve(Amat2,rhs);
detval_dmua_node=rbfemgetdet(phi2_node, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   perturb D at a node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.elemprop=ones(size(cfg.node,1),1);
cfg.elemprop(cfg.elem(perturbeid,1))=2;
cfg.prop(3,:)=cfg.prop(2,:);
cfg.prop(3,2)=cfg.prop(3,2)*1.01;

dD=(1/3/cfg.prop(3,2)-1/3/cfg.prop(2,2));

tic
Amat3=rbfemlhs(cfg); % use native matlab code, 1 sec for 50k nodes
fprintf('solving forward with delta D ... \t%f seconds\n',toc);

phi3_node=rbfemsolve(Amat3,rhs);
detval_dd_node=rbfemgetdet(phi3_node, cfg, loc, bary); % or detval=rbfemgetdet(phi, cfg, rhs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build mua Jacobians, both c/matlab and node/elem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sd=rbsdmap(cfg);

% build all Jacobians using mex code

tic
[Jmua_node_c, Jd_node_c]=rbfemmatrix(cfg, sd, phi);
[Jmua_elem_c, Jd_elem_c]=rbfemmatrix(cfg, sd, phi, cfg.deldotdel, 0);
toc

% build all Jacobians using matlab native code
tic
[Jmua_node_m, Jmua_elem_m, Jd_node_m, Jd_elem_m]=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol); 
fprintf('building J_mua ... \t%f seconds\n',toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Compare mua Jacobian with direct measurement change - c 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dphi_dmua=(detval_dmua-detval)/dmua;          % change of measurement from separate forward
dphi_dmua_node=(detval_dmua_node-detval)/dmua;          % change of measurement from separate forward

dphi_mua_node_c=Jmua_node_c(:,cfg.elem(perturbeid,1));
dphi_mua_node_c=reshape(dphi_mua_node_c,size(cfg.detpos,1),size(cfg.srcpos,1));

dphi_mua_elem_c=Jmua_elem_c(:,perturbeid);  % change of measurement predicted from Jacobians
dphi_mua_elem_c=reshape(dphi_mua_elem_c,size(cfg.detpos,1),size(cfg.srcpos,1));

figure;
len=length(dphi_dmua);
subplot(211);
plot(1:len,dphi_dmua,'r-o',1:len,dphi_mua_elem_c,'b-+');
legend('direct measurement change','predicted from elem J_{\mua} c');
subplot(212);
plot(1:len,dphi_dmua_node,'r-o',1:len,dphi_mua_node_c,'b-+');
legend('direct measurement change','predicted from node J_{\mua} c');

dd_mua_node_c=dphi_mua_node_c./dphi_dmua_node;
fprintf(1,'node-based c Jmua: sum=\t%f\tratio=\t%f\n',sum(Jmua_node_c(:)), median(dd_mua_node_c(:)));
dd_mua_elem_c=dphi_mua_elem_c./dphi_dmua;
fprintf(1,'elem-based c Jmua: sum=\t%f\tratio=\t%f\n',sum(Jmua_elem_c(:)), median(dd_mua_elem_c(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Compare mua Jacobian with direct measurement change - matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dphi_mua_node_m=Jmua_node_m(:,cfg.elem(perturbeid,1));
dphi_mua_node_m=reshape(dphi_mua_node_m,size(cfg.detpos,1),size(cfg.srcpos,1));

dphi_mua_elem_m=Jmua_elem_m(:,perturbeid);  % change of measurement predicted from Jacobians
dphi_mua_elem_m=reshape(dphi_mua_elem_m,size(cfg.detpos,1),size(cfg.srcpos,1));

figure;
len=length(dphi_dmua);
subplot(211);
plot(1:len,dphi_dmua,'r-o',1:len,dphi_mua_elem_m,'b-+');
legend('direct measurement change','predicted from elem J_{\mua} m');
subplot(212);
plot(1:len,dphi_dmua_node,'r-o',1:len,dphi_mua_node_m,'b-+');
legend('direct measurement change','predicted from node J_{\mua} m');

dd_mua_node_m=dphi_mua_node_m./dphi_dmua_node;
fprintf(1,'node-based m Jmua: sum=\t%f\tratio=\t%f\n',sum(Jmua_node_m(:)), median(dd_mua_node_m(:)));
dd_mua_elem_m=dphi_mua_elem_m./dphi_dmua;
fprintf(1,'elem-based m Jmua: sum=\t%f\tratio=\t%f\n',sum(Jmua_elem_m(:)), median(dd_mua_elem_m(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Compare D Jacobian with direct measurement change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dphi_dd=(detval_dd-detval)/dD;          % change of measurement from separate forward
dphi_dd_node=(detval_dd_node-detval)/dD;          % change of measurement from separate forward

%idx=ismember(cfg.elem, cfg.elem(perturbeid,:));
%[ix,iy]=find(idx);
%totalgrad=sum(sum(deldotdel(ix,[1:end 2:4,6:7,9])));
%dphi4=sum(Jd(:,cfg.elem(perturbeid,:)),2)*abs(sum(deldotdel(perturbeid,[1:end 2:4,6:7,9]))/totalgrad);

dphi_d_node_c=Jd_node_c(:,cfg.elem(perturbeid,1));
dphi_d_node_c=reshape(dphi_d_node_c,size(cfg.detpos,1),size(cfg.srcpos,1));

dphi_d_elem_c=Jd_elem_c(:,perturbeid);  % change of measurement predicted from Jacobians
dphi_d_elem_c=reshape(dphi_d_elem_c,size(cfg.detpos,1),size(cfg.srcpos,1));

figure;
subplot(211);
plot(1:len,dphi_dd,'r-o',1:len,dphi_d_elem_c,'b-+');
legend('direct measurement change','predicted from elem J_{D} c');
subplot(212);
plot(1:len,dphi_dd_node,'r-o',1:len,dphi_d_node_c,'b-+');
legend('direct measurement change','predicted from node J_{D} c');

dd_d_node_c=dphi_d_node_c./dphi_dd_node;
fprintf(1,'node-based c Jd: sum=\t%f\tratio=\t%f\n',sum(Jd_node_c(:)), median(dd_d_node_c(:)));
dd_d_elem_c=dphi_d_elem_c./dphi_dd;
fprintf(1,'elem-based c Jd: sum=\t%f\tratio=\t%f\n',sum(Jd_elem_c(:)), median(dd_d_elem_c(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Compare D Jacobian with direct measurement change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dphi_d_node_m=Jd_node_m(:,cfg.elem(perturbeid,1));
dphi_d_node_m=reshape(dphi_d_node_m,size(cfg.detpos,1),size(cfg.srcpos,1));

dphi_d_elem_m=Jd_elem_m(:,perturbeid);  % change of measurement predicted from Jacobians
dphi_d_elem_m=reshape(dphi_d_elem_m,size(cfg.detpos,1),size(cfg.srcpos,1));

figure;
subplot(211);
plot(1:len,dphi_dd,'r-o',1:len,dphi_d_elem_m,'b-+');
legend('direct measurement change','predicted from elem J_{D} m');
subplot(212);
plot(1:len,dphi_dd_node,'r-o',1:len,dphi_d_node_m,'b-+');
legend('direct measurement change','predicted from node J_{D} m');

dd_d_node_m=dphi_d_node_m./dphi_dd_node;
fprintf(1,'node-based m Jd: sum=\t%f\tratio=\t%f\n',sum(Jd_node_m(:)), median(dd_d_node_m(:)));
dd_d_elem_m=dphi_d_elem_m./dphi_dd;
fprintf(1,'elem-based m Jd: sum=\t%f\tratio=\t%f\n',sum(Jd_elem_m(:)), median(dd_d_elem_m(:)));
