function [pointsrc,pointdet,widesrc,widedet]=rbgetoptodes(cfg,wv)
%
% [pointsrc, widesrc]=rbgetoptodes(cfg)
%
% Return the combined list of point optodes (sources+dectors) and widefield
% optodes; In a simulation, all sources, or indenpently, all dectors, can
% only be either point sources or widefield sources.
%
% Edit 2021/11/15 : Now supports simultaneous point source and widefield detectors
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

pointsrc=[];        % EXu - Separated source and detector so easier to concatenate in rbfemrhs
pointdet = [];
widesrc=[];
widedet = [];

if nargin>1
    ltr = rbgetltr(cfg,wv);
else
    ltr = rbgetltr(cfg);
end

if(isfield(cfg,'srcpos') && ~isempty(cfg.srcpos))
    if (size(cfg.srcdir, 1) == size(cfg.srcpos, 1))
        pointsrc=cfg.srcpos+(cfg.srcdir.*ltr);
    elseif (size(cfg.srcdir, 1) == 1)
        pointsrc=cfg.srcpos+repmat(cfg.srcdir.*ltr,size(cfg.srcpos,1),1);
    else
        error('Please provide either one srcdir for all srcpos or one srcdir for each srcpos');
    end
end
if (isfield(cfg,'widesrc') && ~isempty(cfg.widesrc))
    widesrc=cfg.widesrc;
end

if(isfield(cfg,'detpos') && ~isempty(cfg.detpos))
%     pointsrc=[pointsrc; cfg.detpos+repmat(cfg.detdir*ltr,size(cfg.detpos,1),1)];
%     pointdet = cfg.detpos+repmat(cfg.detdir.*ltr,size(cfg.detpos,1),1);
    if (size(cfg.detdir, 1) == size(cfg.detpos, 1))
        pointdet=cfg.detpos+(cfg.detdir.*ltr);
    elseif (size(cfg.detdir, 1) == 1)
        pointdet=cfg.detpos+repmat(cfg.detdir.*ltr,size(cfg.detpos,1),1);
    else
        error('Please provide either one srcdir for all srcpos or one srcdir for each srcpos');
    end
end
if(isfield(cfg,'widedet') && ~isempty(cfg.widedet))
    widedet = cfg.widedet;
end