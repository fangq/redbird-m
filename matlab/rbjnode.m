function [Jmua_node, Jd_node]=rbjnode(Jmua_elem, Jd_elem, elem, nodelen)
if(nargout>nargin)
    error('output number can not exceed input number')
end

Jmua_node=rbelem2node(elem,Jmua_elem,nodelen);

if(nargout>1 && nargin>1)
    Jd_node=rbelem2node(elem,Jd_elem,nodelen);
end