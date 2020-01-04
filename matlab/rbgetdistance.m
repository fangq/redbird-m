function dist=rbgetdistance(srcpos,detpos)
%
% dist=rbgetdistance(srcpos,detpos)
%
% Return the straightline distances between any pair of source and
% detector
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     srcpos: the source position (or boundary condition for widefield src)
%     detpos: the detector position (or boundary condition for widefield det)
%
% output:
%     dist: a Ns x Nd matrix, where Ns is the source number and Nd is the
%           detector number. If either the source or detector is widefield,
%           the dist entries are 0s
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

srcnum=size(srcpos,1);
detnum=size(detpos,1);

if(size(srcpos,2)>4 || size(detpos,2)>4)
    dist=zeros(srcnum,detnum);
    return;
end

dist=repmat(srcpos,detnum,1)-kron(detpos,ones(srcnum,1));
dist=dist.*dist;
dist=sqrt(sum(dist,2));
dist=reshape(dist,srcnum,detnum);