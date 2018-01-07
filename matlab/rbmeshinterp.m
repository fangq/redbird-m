function newval=rbmeshinterp(val,rnode,relem,relemid,relembary)

% interpolate values from one mesh to another
% author: Qianqian Fang, <q.fang at neu.edu>

if(size(val,1)==1)
    val=val(:);
end

newval=zeros(size(rnode,1),size(val,2));
idx=find(~isnan(relemid));
for j=1:size(relembary,2)
    newval(idx,:)=newval(idx,:)+val(relem(relemid(idx),j),:).*repmat(relembary(idx,j),1,size(val,2));
end
