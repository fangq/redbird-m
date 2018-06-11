function nodeval=rbelem2node(elem,elemval,nodelen)
if(isstruct(elem))
    nodelen=size(elem.node,1);
    elem=elem.elem;
end

nelem=size(elem);
nval=size(elemval);

if(nval(2)==nelem(1))
    nodeval=zeros(nval(1),nodelen);
    for j=1:nelem(2)
        nodeval(:,elem(:,j))=nodeval(:,elem(:,j))+elemval(:,elem(:,j));
    end
else
    nodeval=zeros(nodelen,nval(1));
    for j=1:nelem(2)
        nodeval(elem(:,j),:)=nodeval(elem(:,j),:)+elemval(elem(:,j),:);
    end
end
nodeval=nodeval*0.25;