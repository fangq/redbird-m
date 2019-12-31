function optode=rbgetoptodes(cfg)

optode=[];

if(isfield(cfg,'srcpos') && ~isempty(cfg.srcpos))
    if(size(cfg.srcpos,2) ~= size(cfg.face,1))
	optode=cfg.srcpos(:,1:3);
    end
end

if(isfield(cfg,'detpos') && ~isempty(cfg.detpos))
    optode=[optode; cfg.detpos(:,1:3)];
end