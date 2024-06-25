function dist=rbgetdistance(srcpos,detpos,badsrc,baddet,widesrc,widedet,cfg)
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

srcnum=size(srcpos,1) + size(widesrc,1);
detnum=size(detpos,1) + size(widedet,1);

if(size(srcpos,2)>4 || size(detpos,2)>4)
    dist=zeros(srcnum,detnum);
    return;
end

wsrcpos = [];
wdetpos = [];

if ~isempty(widesrc)
    wsrcpos = rbcomputewfcom(widesrc,cfg);
    srcpos = [srcpos; wsrcpos];
end
if ~isempty(widedet)
    wdetpos = rbcomputewfcom(widedet,cfg);
    detpos = [detpos; wdetpos];
end

if(nargin<4)
    baddet=[];
    if(nargin<3)
        badsrc=[];
    end
end

goodsrc=setdiff(1:srcnum,badsrc);
gooddet=setdiff(1:detnum,baddet);

dd=repmat(srcpos(goodsrc,:),length(gooddet),1)-kron(detpos(gooddet,:),ones(length(goodsrc),1));
dd=dd.*dd;
dd=sqrt(sum(dd,2));
dd=reshape(dd,length(goodsrc),length(gooddet));

dist=inf*ones(srcnum,detnum);
dist(goodsrc,gooddet)=dd;