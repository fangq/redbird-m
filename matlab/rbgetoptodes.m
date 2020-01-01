function [pointsrc, widesrc]=rbgetoptodes(cfg)

pointsrc=[];
widesrc=[];

if(isfield(cfg,'srcpos') && ~isempty(cfg.srcpos))
    if(size(cfg.srcpos,2) == size(cfg.face,1))
        widesrc=cfg.srcpos;
    else
        pointsrc=cfg.srcpos;
    end
end

if(isfield(cfg,'detpos') && ~isempty(cfg.detpos))
    if(size(cfg.detpos,2) == size(cfg.face,1))
        widesrc=[widesrc; cfg.detpos];
    else
        pointsrc=[pointsrc; cfg.detpos];
    end
end