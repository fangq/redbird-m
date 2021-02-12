function Lmat=rbprior(seg,priortype,param)

if(nargin<2)
    priortype='';
end

if(isempty(priortype))
    Lmat=[];
    return;
end

if(isvector(seg))
    if(isempty(priortype) && all(mod(seg(:),1) == 0))
        [labels,ix,iy]=unique(seg);
        clear ix;
        counts=hist(seg,labels);
        Lmat=spzeros(length(seg));
        if(strcmp(priortype,'laplace'))
            for i=1:length(labels)
                idx=find(iy==labels(i));
                if(counts>0)
                    Lmat(idx,idx)=-1/counts(i);
                end
            end
            Lmat(1:size(Lmat,1)+1:end)=1;
            return;
        elseif(strcmp(priortype,'helmholtz') && nargin>2)
            for i=1:length(labels)
                idx=find(iy==labels(i));
                if(counts>0)
                    Lmat(idx,idx)=-1./(counts(i)+param(idx,idx));
                end
            end
            Lmat(1:size(Lmat,1)+1:end)=1;
            return;
        end
    end
    if(all((seg<=1.0) && (seg>=0)))
        seg=[seg(:) 1-seg(:)];
    end
end

if(all((seg(:)<=1.0) & (seg(:)>=0)) && size(seg,2)>1 && strcmp(priortype,'comp'))
    spidx=sparse(1:length(seg));
    cdist=@(a,b) compdist(a,b,seg,param);
    Lmat=bsxfun(cdist,spidx,spidx');
    rowsum=abs(sum(Lmat)).';
    [ix,iy]=find(Lmat);
    Lmat(find(Lmat))=Lmat(find(Lmat))./(param.beta*sqrt(rowsum(ix).*rowsum(iy)));
    Lmat=Lmat+speye(length(seg));
    return;
end

function val=compdist(a,b,seg,param)
dval=sum(abs(repmat(seg(a,:),size(b,1),1)-seg(b,:)),2);
dval=dval.*dval;
val=sparse((dval<(param.alpha*size(seg,2))*(param.alpha*size(seg,2))).*dval);
if(any(val>0))
    val(val>0)=-param.alpha-sqrt(dval(val>0))/size(seg,2);
end
