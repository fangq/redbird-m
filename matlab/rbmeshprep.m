function newcfg=rbmeshprep(cfg)

if(~isfield(cfg,'node') || ~isfield(cfg,'elem'))
    error('cfg.node or cfg.elem is missing');
end
if(~isfield(cfg,'elemprop') ||isempty(cfg.elemprop) && size(cfg.elem,2)>4)
    cfg.elemprop=cfg.elem(:,5);
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

if(~isfield(cfg,'cols') || isempty(cfg.cols))
    [cfg.rows,cfg.cols,cfg.idxcount]=rbfemnz(cfg.elem,size(cfg.node,1));
end

if(~isfield(cfg,'idxsum') || isempty(cfg.idxsum))
    cfg.idxsum=cumsum(cfg.idxcount);
end

newcfg=cfg;