function [Jd_node, Jd_elem]=rbjacdcoef(sd, phi, deldotdel, felem)

if(nargin<4 || isempty(sd) || isempty(phi) || isempty(deldotdel))
    error('you must give at least 5 inputs and they must not be empty');
end

Jd_node=zeros(size(sd,1),size(phi,1));

nelem=size(deldotdel,1);
Jd_elem=zeros(size(sd,1),nelem);

idx=[1 1 1 1 2 2 2 3 3 4; 1 2 3 4 2 3 4 3 4 4];
idx2=[2 3 4 3 4 4; 1 1 1 2 2 3; 2 3 4 6 7 9];

for i=1:nelem
    Jcol=sum(phi(felem(i,idx(1,:)),sd(:,1)).*phi(felem(i,idx(2,:)),sd(:,2)).*repmat(deldotdel(i,:)',1,size(sd,1)),1).';
    Jcol=Jcol+sum(phi(felem(i,idx2(1,:)),sd(:,1)).*phi(felem(i,idx2(2,:)),sd(:,2)).*repmat(deldotdel(i,idx2(3,:))',1,size(sd,1)),1).';
    Jd_node(:,felem(i,:))=Jd_node(:,felem(i,:))+repmat(Jcol,1,4);
    Jd_elem(:,i)=Jcol;
end
Jd_node=-Jd_node*0.25;
Jd_elem=-Jd_elem;
