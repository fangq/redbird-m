function ltr=rbgetltr(cfg,wavelength)

if(isfield(cfg,'bulk') && isfield(cfg.bulk,'musp'))
    if(isfield(cfg.bulk,'mua'))
        ltr=1/(3*(cfg.bulk.musp+cfg.bulk.mua));
    else
        ltr=1/(3*cfg.bulk.musp);
    end
else
    if(isa(cfg.prop,'containers.Map'))
        if(nargin>1)
            prop=cfg.prop(wavelength); % should 
        else
            wv=cfg.prop.keys;
            prop=cfg.prop(wv{1}); % should
        end
    else
        prop=cfg.prop;
    end
    if(size(prop,2)>2)
        ltr=1/(prop(2,1)+prop(2,2)*(1-prop(2,3)));
    else
        ltr=1/(prop(2,1)+prop(2,2));
    end
end
