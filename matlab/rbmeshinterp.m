function newval=rbmeshinterp(val,relemid,relembary,fromelem)

% interpolate values from one mesh to another
% author: Qianqian Fang, <q.fang at neu.edu>

if(size(val,1)==1)
    val=val(:);
end

newval=zeros(length(relemid),size(val,2));

idx=find(~isnan(relemid));

for i=1:size(val,2)
    newval(idx,:)=newval(idx,:) + reshape(val(fromelem(relemid(idx),:),i),1,size(relembary,2)).*relembary(idx,:);
end