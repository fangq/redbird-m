function optode=rbgetoptodes(cfg)

optode=[];

if(isfield(cfg,'srcpos') && ~isempty(cfg.srcpos))
    optode=cfg.srcpos(:,1:3);
end

if(isfield(cfg,'detpos') && ~isempty(cfg.detpos))
    optode=[optode; cfg.detpos(:,1:3)];
end