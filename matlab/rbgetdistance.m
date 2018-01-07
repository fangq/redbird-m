function dist=rbgetdistance(srcpos,detpos)

srcnum=size(srcpos,1);
detnum=size(detpos,1);

dist=repmat(srcpos,detnum,1)-kron(detpos,ones(srcnum,1));
dist=dist.*dist;
dist=sqrt(sum(dist,2));
dist=reshape(dist,srcnum,detnum);