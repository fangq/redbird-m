function sd=rbsdmap(cfg,maxdist)

if(~isfield(cfg,'srcpos') || isempty(cfg.srcpos) || ~isfield(cfg,'detpos') || isempty(cfg.detpos))
    error('you must define at least 1 source and 1 detector in order to use this function');
end

dist=rbgetdistance(cfg.srcpos,cfg.detpos);

srcnum=size(cfg.srcpos,1);
detnum=size(cfg.detpos,1);

[ss,dd]=meshgrid(1:srcnum,srcnum+1:srcnum+detnum);
sd=[ss(:),dd(:)];
if(nargin<2)
    sd(:,3)=1;
else
    sd(:,3)=dist(:)<=maxdist;
end