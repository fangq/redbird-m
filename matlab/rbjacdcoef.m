function [Jd, JD]=rbjacdcoef(sd, phi, deldotdel, felem)

if(nargin<4 || isempty(sd) || isempty(phi) || isempty(deldotdel))
    error('you must give at least 5 inputs and they must not be empty');
end

Jd=zeros(size(sd,1),size(phi,1));

idx=[2 3 4 6 7 9
     1 1 1 2 2 3
     2 3 4 3 4 4];
id0=[1 5 8 10];
 
nmax=max(felem);
nelem=size(deldotdel,1);

JD=zeros(size(sd,1),nelem);

for i=1:size(sd,1)
    JD(i,:)=JD(i,:)+sum(deldotdel(:,idx(1,:)).*reshape(phi(felem(:,idx(2,:)),sd(i,1)).*phi(felem(:,idx(3,:)),sd(i,2))+...
                    phi(felem(:,idx(2,:)),sd(i,2)).*phi(felem(:,idx(3,:)),sd(i,1)),nelem,size(idx,2)),2)';
    for j=1:size(id0,2)
        JD(i,:)=JD(i,:)+(deldotdel(:,id0(j)).*phi(felem(:,j),sd(i,1)).*phi(felem(:,j),sd(i,2)))';
    end
    for j=1:4
        val=accumarray(felem(:,j),JD(i,:))';
        Jd(i,1:nmax(j))=Jd(i,1:nmax(j))+val;
    end
end
Jd=-Jd*0.25;
JD=-JD;
