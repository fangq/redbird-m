function newval=rbmeshinterp(fromval,elemid,elembary,fromelem)

% interpolate values from one mesh to another
% author: Qianqian Fang, <q.fang at neu.edu>

if(size(fromval,1)==1)
    fromval=fromval(:);
end

idx=find(~isnan(elemid));

allval=reshape(fromval(fromelem(elemid(idx),:),:),length(idx),size(elembary,2),size(fromval,2));
tmp=cellfun(@(x) sum(elembary(idx,:).*x,2), num2cell(allval,[1 2]),'UniformOutput',false);
newval=cat(3,tmp{:});

