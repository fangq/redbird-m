function cfg=rbsetmesh(cfg0,node,elem)

names=fieldnames(cfg0);
names=intersect(names,{'face','evol','deldotdel','isreoriented','nvol','cols','idxsum','elemprop'});

cfg0.node=node;
cfg0.elem=elem;

if(~isempty(names))
    cfg=rmfield(cfg0,names);
else
    cfg=cfg0;
end
cfg=rbmeshprep(cfg);
