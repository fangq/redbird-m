function prop=rbupdateprop(cfg)
%
% prop=rbupdateprop(cfg)
%
% Update the direct material properties (cfg.prop - optical or EM 
% properties used by the forward solver) using drived properties 
% (cfg.param - physiological parameters) at give wavelengths (cfg.wavelen)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird data structure 
%
% output:
%     prop: the updated material property data to replace cfg.prop
%
% license:
%     BSD or GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(~isfield(cfg,'prop') || ~isa(cfg.prop,'containers.Map') ||  ~isfield(cfg,'param')))
    error('input cfg must be a struct and must have subfield names "prop" and "param"');
end

if(~isa(cfg.prop,'containers.Map'))
    prop=cfg.prop;
    return;
end

wv=keys(cfg.prop);

for i=1:length(wv)
    wavelen=wv{i};
    types=fieldnames(cfg.param);
    types=intersect(cfg.param,{'hbo','hbr','water','lipids','aa3'});
    if(isempty(types))
        error('specified parameters are not supported');
    end
    extin=rbextinction(str2num(wv{i}), types);
    mua=zeros(size(cfg.param(types{1})));
    for j=1:length(types)
        mua=mua+extin.(types{i})*cfg,param(types{i});
    end
    musp=cfg.param('scatamp').*exp(-cfg.wavelen(i).*cfg.param('scatpower'));
    prop(wavelen)=[mua musp];
end
