function newdata=rbmasksum(data, mask)
dim=size(mask);
if(dim(1)==1 && dim(2)>1 && length(mask)==size(data,2))
    [xi,yi]=ndgrid(1:size(data,1),mask(:));
else
    [xi,yi]=ndgrid(mask(:),1:size(data,2));
end
newdata=accumarray([xi(:),yi(:)],data(:));