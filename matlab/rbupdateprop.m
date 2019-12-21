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
    musp=cfg.param('scatamp').*exp(-cfg.wavelen(i).*cfg.param('scatpower'));
    prop.(wavelen)=[mua musp];
end