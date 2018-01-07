function [Amat,deldotdel]=rbfemlhs(cfg)

% create the FEM stiffness matrix (left-hand-side) for solving the diffusion equation

nn=size(cfg.node,1);
[Adiag, Aoffd, deldotdel]=rbfemmatrix(cfg);
Amat = sparse([cfg.rows,cfg.cols,(1:nn)],[cfg.cols,cfg.rows,(1:nn)],[Aoffd,Aoffd,Adiag],nn,nn);
deldotdel=deldotdel';
