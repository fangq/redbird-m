function dist=rbgetdistance(srcpos,detpos)

srcnum=size(srcpos,1);
detnum=size(detpos,1);

if(size(srcpos,2)>4 || size(detpos,2)>4)
    dist=zeros(srcnum,detnum);
    return;
end

dist=repmat(srcpos,detnum,1)-kron(detpos,ones(srcnum,1));
dist=dist.*dist;
dist=sqrt(sum(dist,2));
dist=reshape(dist,srcnum,detnum);