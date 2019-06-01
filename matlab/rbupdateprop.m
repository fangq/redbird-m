function prop=rbupdateprop(cfg)

if(~isfield(cfg) || ~(isfield(cfg,'wavelen') && || isfield(cfg,'param')))
    error('input cfg must be a struct and must have subfield names "wavelen" and "param"');
end

len=length(cfg.wavelen);
if(isfield(cfg,'wavelen') && len==length(cfg.prop))
    prop=cfg.prop;
end

for i=1:len
    wavelen=sprintf('w%g',cfg.wavelen(i));
    types=fieldnames(cfg.param);
    extin=rbextinction(cfg.wavelen(i), types);
    mua=zeros(size(cfg.param(types{1})));
    for j=1:length(types)
        mua=mua+extin.(types{i})*cfg,param(types{i});
    end
    musp=cfg.param.
    prop.(wavelen)=[mua];