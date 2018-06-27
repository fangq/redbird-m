function newval=rbmeshremap(fromval,elemid,elembary,toelem,nodeto)
%
% newval=rbmeshremap(fromval,elemid,elembary,toelem,nodeto)
%
% Redistribute nodal values from the source mesh to the target mesh so that 
% the sum of each property on each mesh is the same
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%	 fromval: values defined at the source mesh nodes, the row or column
%	          number must be the same as the source mesh node number, which
%	          is the same as the elemid length
%	 elemid: the IDs of the target mesh element that encloses the nodes of
%            the source mesh nodes; a vector of length of src mesh node
%            count; elemid and elembary can be generated by calling
%
%           [elemid,elembary]=tsearchn(node_target, elem_target, node_src);
%
%           note that the mapping here is inverse to that in meshinterp()
%
%	 elembary: the bary-centric coordinates of each source mesh nodes
%	         within the target mesh elements, sum of each row is 1, expect
%	         3 or 4 columns (or can be N-D)
%    toelem: the element list of the target mesh
%    nodeto: the total number of target mesh nodes
%
%
% output:
%	 newval: a 2D array with rows equal to the target mesh nodes (nodeto), 
%            and columns equals to the value numbers defined at each source
%            mesh node
% example:
%
%    [n1,f1,e1]=meshabox([0 0 0],[10 20 5],1); % src mesh
%    [n2,f2,e2]=meshabox([0 0 0],[10 20 5],2); % target mesh
%    [id, ww]=tsearchn(n2,e2,n1);              % project src to target mesh
%    value_src=n1(:,[2 3 1]);             % create dummy values at src mesh
%    newval=rbmeshremap(value_src,id,ww,e2,size(n2,1)); % map to target
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(size(fromval,1)==1)
    fromval=fromval(:);
end

if(size(fromval,2)==length(elemid))
    fromval=fromval.';
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