function newval=rbmeshremap(fromval,elemid,elembary,toelem,nodeto)

% redistribute values from one mesh to another so that the sum is the same
% author: Qianqian Fang, <q.fang at neu.edu>

if(size(fromval,1)==1)
    fromval=fromval(:);
end

if(size(formval,2)==length(elemid))
    formval=formval.';
end

newval=zeros(nodeto,size(fromval,2));

idx=~isnan(elemid);
idx=elemid(idx);

nodeval=repmat(fromval,1,1,size(elembary,2)).*repmat(permute(elembary,[1,3,2]),1,size(fromval,2),1);

for i=1:size(elembary,2)
    [ix,iy]=meshgrid(toelem(idx,i),1:size(fromval,2));
    nval=nodeval(:,:,i).';
    newval=newval + accumarray([ix(:),iy(:)],nval(:), size(newval));
end