function [deldotdel, delphi]=rbdeldotdel(cfg)
% Format:
%   [deldotdel, delphi]=rbdeldotdel(cfg)
%
% Compute deldotdel=<grad(phi_i).grad(phi_j)>, where <> means spatial integration
% inside elements, "." means dot-product, grad() means gradience, phi means
% linear basis function in a tetrahedron. For a linear function phi,
% grad(phi) is a constant across the element.
%
% Author:
%   Qianqian Fang <q.fang at neu.edu>
%

no=cfg.node;
el=cfg.elem(:,1:4);

no=reshape(no(el',:)',3,4,size(el,1));

delphi=zeros(3,4,size(el,1));

col=[4 2 3 2
     3 1 4 3
     2 4 1 4
     1 3 2 1];

for coord=1:3
    idx=1:3;
    idx(coord)=[];
    for i=1:4
        delphi(coord,i,:)=squeeze(((no(idx(1),col(i,1),:)-no(idx(1),col(i,2),:)).*(no(idx(2),col(i,3),:)-no(idx(2),col(i,4),:))-...
                                   (no(idx(1),col(i,3),:)-no(idx(1),col(i,4),:)).*(no(idx(2),col(i,1),:)-no(idx(2),col(i,2),:))))./(cfg.evol(:)*6);
    end
end

deldotdel=zeros(size(el,1),10);
count=1;

for i=1:4
    for j=i:4
        deldotdel(:,count)=sum(squeeze(delphi(:,i,:).*delphi(:,j,:)),1);
        count=count+1;
    end
end
deldotdel=deldotdel.*repmat(cfg.evol(:),1,10);