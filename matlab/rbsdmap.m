function [sd,dist]=rbsdmap(cfg,maxdist)
%
% [sd,dist]=rbsdmap(cfg,maxdist)
%
% Create a source-detector mapping table (sd) by setting a maximum separation distance

% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird data structure
%     maxdist: the maximum source detector separation distance (in mm)
%
% output:
%     sd: the source detector mapping table
%     dist: the full table of all source detector separations (#row=#sources, #col=#detectors)
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

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