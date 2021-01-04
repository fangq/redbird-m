function prop=rbupdateprop(cfg,wv)
%
% prop=rbupdateprop(cfg)
%
% Update the direct material properties (cfg.prop - optical or EM 
% properties used by the forward solver) using derived properties 
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

% keys(cfg.param) provides chromorphore species, 2nd input wv or 
% keys(cfg.prop) provides the wavelength list, thus, both must be
% containers.Map

if(~isfield(cfg,'prop') || ~isa(cfg.prop,'containers.Map') ||  ~isfield(cfg,'param'))
    error('input cfg must be a struct and must have subfield names "prop" and "param"');
end

if(nargin<2)
    wv=keys(cfg.prop);
end

prop=containers.Map();

for i=1:length(wv)
    wavelen=wv{i};
    types=intersect(fieldnames(cfg.param),{'hbo','hbr','water','lipids','aa3'});
    if(isempty(types))
        error('specified parameters are not supported');
    end
    extin=rbextinction(str2double(wv{i}), types);
    mua=zeros(size(cfg.param.(types{1})));
    for j=1:length(types)
        mua=mua+extin(i)*cfg.param.(types{i});
    end
    if(isfield(cfg.param,'scatamp') && isfield(cfg.param,'scatpower'))
        musp=cfg.param.('scatamp').*exp(-cfg.wavelen(i).*cfg.param.('scatpower'));
        prop(wavelen)=[mua musp];
    else
        prop(wavelen)=mua;
    end 
end
