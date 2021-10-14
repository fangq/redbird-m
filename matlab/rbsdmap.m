function [sd,dist]=rbsdmap(cfg,varargin)
%
% [sd,dist]=rbsdmap(cfg,maxdist)
% [sd,dist]=rbsdmap(cfg,'maxdist',10,'excludesrc',[1,2],'excludedet',[1],'wavesrc',...,'wavedet',...)
%
% Create a source-detector mapping table (sd) by setting a maximum separation distance
%
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

if(length(varargin)==1)
    maxdist=varargin{1};
elseif(isempty(varargin))
    maxdist=inf;
elseif(length(varargin)>1)
    opt=varargin2struct(varargin{:});
    maxdist=jsonopt('maxdist',inf,opt);
end

if ~exist('opt','var')
    opt = struct();
end

srcnum=size(cfg.srcpos,1);
detnum=size(cfg.detpos,1);

badsrc=jsonopt('excludesrc',[],opt);
baddet=jsonopt('excludedet',[],opt);

if(isfield(cfg,'srcpos') && ((size(cfg.srcpos,2) == size(cfg.face,1)) || (size(cfg.srcpos,2) == size(cfg.node,1))) )
    dist=zeros(srcnum,detnum);
else
    dist=rbgetdistance(cfg.srcpos,cfg.detpos,badsrc,baddet);
end

goodsrc=sort(setdiff(1:srcnum,badsrc));
gooddet=sort(setdiff(1:detnum,baddet));

if(isfield(cfg,'prop') && isa(cfg.prop,'containers.Map'))
    wavelengths=cfg.prop.keys;
    sd=containers.Map();
    [ss,dd]=meshgrid(goodsrc,srcnum+gooddet);

    for wv=wavelengths
        wid=wv{1};
        if((isfield(opt,'wavesrc') && ~isempty(opt.wavesrc)) || ...
           (isfield(opt,'wavedet') && ~isempty(opt.wavedet)))
             wavesrc=jsonopt('wavesrc',containers.Map({wid},{[]}), opt);
             if(~isempty(wavesrc(wid)))
                wavesrc = intersect(goodsrc,wavesrc(wid));
             else
                wavesrc=setdiff(goodsrc,wavesrc(wid));
             end
             wavedet=jsonopt('wavedet',containers.Map({wid},{[]}), opt);
             if(~isempty(wavedet(wid)))
                wavedet = intersect(goodsrc,wavedet(wid));
             else
                wavedet=setdiff(goodsrc,wavedet(wid));
             end
             [ss,dd]=meshgrid(wavesrc,srcnum+wavedet);
        elseif ((isfield(cfg,'wavesrc') && ~isempty(cfg.wavesrc)) || ...
           (isfield(cfg,'wavedet') && ~isempty(cfg.wavedet)))
             if(~isempty(cfg.wavesrc(wid)))
                wavesrc = intersect(goodsrc,cfg.wavesrc(wid));
             else
                wavesrc = setdiff(goodsrc,cfg.wavesrc(wid));
             end
             if(~isempty(cfg.wavedet(wid)))
                wavedet = intersect(goodsrc,cfg.wavedet(wid));
             else
                wavedet = setdiff(goodsrc,cfg.wavedet(wid));
             end
             [ss,dd]=meshgrid(wavesrc,srcnum+wavedet);
        end
        sdwv=[ss(:),dd(:)];
%         if(nargin<2 || (size(cfg.srcpos,2) == size(cfg.face,1)))
        if( isinf(maxdist) )
            sdwv(:,3)=1;
        else
            sdwv(:,3)=(dist(:)<maxdist);
        end
        
        if (isfield(cfg,'rfcw'))
            modes = cfg.rfcw.src.keys;
            for md = modes
                mid = md{1};
                mdSRC = cfg.rfcw.src(mid);
                mdDET = cfg.rfcw.det(mid) + srcnum;
                
                mdChan = (ismember(sdwv(:,1), mdSRC) & ismember(sdwv(:,2),mdDET));
                if strcmp(mid,'RF')
                    sdwv(mdChan,4) = 1;
                elseif strcmp(mid,'CW')
                    sdwv(mdChan,4) = 2;
                end
            end
        end
        sd(wid)=sdwv;
    end
else
    [ss,dd]=meshgrid(goodsrc,srcnum+gooddet);
    sd=[ss(:),dd(:)];
    if(nargin<2 || (size(cfg.srcpos,2) == size(cfg.face,1)))
        sd(:,3)=1;
    else
        sd(:,3)=(dist(:)<maxdist);
    end
end