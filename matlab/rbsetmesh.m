function cfg=rbsetmesh(cfg0,node,elem,prop,propidx)

names=fieldnames(cfg0);
names=intersect(names,{'face','evol','deldotdel','isreoriented','nvol','cols','idxsum','musp0'});

cfg0.node=node;
cfg0.elem=elem;
if(nargin>3)
    cfg0.prop=prop;
end
if(nargin>4)
    cfg0.seg=propidx;
end

if(~isempty(names))
    cfg=rmfield(cfg0,names);
else
    cfg=cfg0;
end
cfg=rbmeshprep(cfg);
