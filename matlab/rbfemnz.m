function [rows,cols,connnum]=rbfemnz(elem,nn)
% rbfemnz: store the indices of all off-diagonal entries of FEM matrix
% author: fangq (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/11/21
%
% parameters:
%    elem:  element table of a mesh
%    nn:    number of all nodes
%    rows: output, the row index of the non-zero off-diagonal elements
%    cols: output, the column index of the non-zero off-diagonal elements
%    connnum: output, total non-zero off-diagonal element number

[conn,connnum,count]=meshconn(elem,nn);

rows=zeros(1,count);
cols=zeros(1,count);
count=0;
for i=1:nn
    rows(count+1:count+connnum(i))=i;
    cols(count+1:count+connnum(i))=conn{i};
    count=count+connnum(i);
end
