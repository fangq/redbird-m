function [Jmua_node, Jmua_elem, Jd_node, Jd_elem]=rbjac(sd, phi, deldotdel, felem, evol)

if(nargin<5 || isempty(sd) || isempty(phi) || isempty(deldotdel)|| isempty(evol))
    error('you must give at least 5 inputs and they must not be empty');
end

nelem=size(felem,1);

Jmua_node=zeros(size(sd,1),size(phi,1));
Jmua_elem=zeros(size(sd,1),nelem);
if(nargout>2)
    Jd_node=zeros(size(sd,1),size(phi,1));
    Jd_elem=zeros(size(sd,1),nelem);
end

idx=[1 1 1 2 2 3 2 3 4 3 4 4
     2 3 4 3 4 4 1 1 1 2 2 3
     2 3 4 6 7 9 2 3 4 6 7 9];

for i=1:nelem
    phidotphi1=phi(felem(i,1:4),sd(:,1)).*phi(felem(i,1:4),sd(:,2));
    phidotphi2=phi(felem(i,idx(1,:)),sd(:,1)).*phi(felem(i,idx(2,:)),sd(:,2));
    Jmua_elem(:,i)=-(sum(phidotphi1,1)+sum(phidotphi2,1)*0.5)*(0.1*evol(i));
    if(nargout>2)
        Jcol=phidotphi1.'*deldotdel(i,[1 5 8 10]).';
        Jcol=Jcol+phidotphi2.'*deldotdel(i,idx(3,:)).';
        Jd_elem(:,i)=-Jcol;
    end
    for j=1:4
        Jmua_node(:,felem(i,j))=Jmua_node(:,felem(i,j))+Jmua_elem(:,i);
    end
    if(nargout>2)
        for j=1:4
            Jd_node(:,felem(i,j))=Jd_node(:,felem(i,j))-Jcol;
        end
    end
end
Jmua_node=Jmua_node*0.25;
if(nargout>2)
    Jd_node=Jd_node*0.25;
end