function newval=rbmeshremap(val,relemid,relembary,toelem,nodeto)

% redistribute values from one mesh to another so that the sum is the same
% author: Qianqian Fang, <q.fang at neu.edu>

if(size(val,1)==1)
    val=val(:);
end

newval=zeros(nodeto,size(val,2));

idx=~isnan(relemid);
idx=relemid(idx);

nodeval=repmat(val,1,1,size(relembary,2)).*repmat(permute(relembary,[1,3,2]),1,size(val,2),1);

for i=1:size(relembary,2)
    [ix,iy]=meshgrid(toelem(idx,i),1:size(val,2));
    nval=nodeval(:,:,i).';
    newval=newval + accumarray([ix(:),iy(:)],nval(:), size(newval));
end