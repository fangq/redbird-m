function newcfg=rbmeshprep(cfg)
%
% newcfg=rbmeshprep(cfg)
%
% Compute all missing fields from the cfg input sturcture to get 
% ready for forward and inverse modeling
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the initial simulation data structure
%
% output:
%     newcfg: the updated simulation data structure after adding all missing fields
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(~isfield(cfg,'node') || ~isfield(cfg,'elem'))
    error('cfg.node or cfg.elem is missing');
end

cfg.elem(:,1:4)=meshreorient(cfg.node(:,1:3),cfg.elem(:,1:4));

if((~isfield(cfg,'seg') ||isempty(cfg.seg)) && size(cfg.elem,2)>4)
    cfg.seg=cfg.elem(:,5);
    cfg.elem(:,5)=[];
end
if(~isfield(cfg,'isreoriented') || isempty(cfg.isreoriented) || cfg.isreoriented==0)
    cfg.elem=meshreorient(cfg.node,cfg.elem(:,1:4));
    cfg.isreoriented=1;
end
if(~isfield(cfg,'face') || isempty(cfg.face))
    cfg.face=volface(cfg.elem);
end
if(~isfield(cfg,'area') || isempty(cfg.area))
    cfg.area=elemvolume(cfg.node,cfg.face);
end
if(~isfield(cfg,'evol') || isempty(cfg.evol))
    cfg.evol=elemvolume(cfg.node,cfg.elem);
end
if(~isfield(cfg,'nvol') || isempty(cfg.nvol))
    cfg.nvol=nodevolume(cfg.node,cfg.elem, cfg.evol);
end
if(find(cfg.evol==0))
    fprintf(1,['degenerated elements are detected: [' sprintf('%d ',find(cfg.evol==0)) ']\n']);
    error(['input mesh can not contain degenerated elements, ' ...
        'please double check your input mesh; if you use a ' ...
        'widefield source, please rerun mmcsrcdomain and setting ' ...
        '''Expansion'' option to a larger value (default is 1)']);
end
if(~isfield(cfg,'srcpos'))
    error('cfg.srcpos field is missing');
end
if(~isfield(cfg,'srcdir'))
    error('cfg.srcdir field is missing');
end
if(isfield(cfg,'prop') && isfield(cfg,'param') && ...
        isa(cfg.prop,'containers.Map'))
    wv=cfg.prop.keys;
    if(~isempty(wv))
        cfg.prop=rbupdateprop(cfg);
    end
end
% compute R_eff - effective reflection coeff, and musp0 - background mus'
if(~isfield(cfg,'reff') || isempty(cfg.reff))
    bkprop=rbgetbulk(cfg);
    if(isa(bkprop,'containers.Map'))
        cfg.reff=containers.Map();
        cfg.musp0=containers.Map();
        for waveid=bkprop.keys
            wv=waveid{1};
            prop=bkprop(wv);
            cfg.reff(wv)=rbgetreff(prop(4), 1);
            cfg.musp0(wv)=prop(2)*(1-prop(3));
        end
    else
        cfg.reff=rbgetreff(bkprop(4), 1);
        cfg.musp0=bkprop(2)*(1-bkprop(3));
    end
end
if(isfield(cfg,'srctype') && ~ismember(cfg.srctype,{'pencil','isotropic'}))
    cfg.srcpos0=cfg.srcpos;
    cfg.srcpos=rbsrc2bc(cfg);
end
if(isfield(cfg,'dettype') && ~ismember(cfg.dettype,{'pencil','isotropic'}))
    cfg.detpos0=cfg.detpos;
    cfg.detpos=rbsrc2bc(cfg,1);
end
if(~isfield(cfg,'cols') || isempty(cfg.cols))
    [cfg.rows,cfg.cols,cfg.idxcount]=rbfemnz(cfg.elem,size(cfg.node,1));
end

if(~isfield(cfg,'idxsum') || isempty(cfg.idxsum))
    cfg.idxsum=cumsum(cfg.idxcount);
end

if(~isfield(cfg,'deldotdel') || isempty(cfg.deldotdel))
    cfg.deldotdel=rbdeldotdel(cfg);
end
newcfg=cfg;