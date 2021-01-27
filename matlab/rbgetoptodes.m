function [pointsrc, widesrc]=rbgetoptodes(cfg)
%
% [pointsrc, widesrc]=rbgetoptodes(cfg)
%
% Return the combined list of point optodes (sources+dectors) and widefield
% optodes; In a simulation, all sources, or indenpently, all dectors, can
% only be either point sources or widefield sources.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the initial simulation data structure
%
% output:
%     pointsrc: combined point source list - a dimension of Np x 3, where
%          Np is the total point source+point detector number
%     widesrc: combined widefield source list - a dimension of Nw x 3, where
%          Nw is the total point source+point detector number
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

pointsrc=[];
widesrc=[];

ltr=rbgetltr(cfg);

if(isfield(cfg,'srcpos') && ~isempty(cfg.srcpos))
    if(size(cfg.srcpos,2) == size(cfg.node,1))
        widesrc=cfg.srcpos;
    else
        pointsrc=cfg.srcpos+repmat(cfg.srcdir*ltr,size(cfg.srcpos,1),1);
    end
end

if(isfield(cfg,'detpos') && ~isempty(cfg.detpos))
    if(size(cfg.detpos,2) == size(cfg.node,1))
        widesrc=[widesrc; cfg.detpos];
    else
        pointsrc=[pointsrc; cfg.detpos+repmat(cfg.detdir*ltr,size(cfg.detpos,1),1)];
    end
end