function [detval, goodidx]=rbfemgetdet(phi, cfg, optodeloc, optodebary)

% output:
%   detval: #det x #src array, denoting the measurement data

if(~isfield(cfg,'srcpos') || isempty(cfg.srcpos) || ~isfield(cfg,'detpos') || isempty(cfg.detpos))
    error('you must define at least 1 source and 1 detector in order to use this function');
end

srcnum=size(cfg.srcpos,1);
detnum=size(cfg.detpos,1);

goodidx=find(~isnan(optodeloc(srcnum+1:srcnum+detnum)));
detval=zeros(detnum,srcnum);
gooddetval=zeros(length(goodidx),srcnum);

if(nargin==3)
    detval=optodeloc(:,srcnum+1:srcnum+detnum)'*phi(:,1:srcnum);
else
    for i=1:length(goodidx)
        if(~isnan(optodeloc(i)))
            gooddetval(i,:)=sum(phi(cfg.elem(optodeloc(srcnum+goodidx(i)),:),1:srcnum).*repmat(optodebary(srcnum+goodidx(i),:)',1,srcnum),1);
        end
    end
    detval(goodidx,:)=gooddetval;
end