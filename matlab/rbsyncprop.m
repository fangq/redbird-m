function [cfg, recon] = rbsyncprop(cfg, recon)

if (isfield(recon, 'node') && isfield(recon, 'elem')) % dual-mesh
    labelmax = min(size(recon.node, 1), size(recon.elem, 1));
else % single-mesh
    labelmax = min(size(cfg.node, 1), size(cfg.elem, 1));
end

if (isfield(recon, 'param')) % map recon.param to cfg.param
    if (~isstruct(recon.param))
        error('recon.param must be a struct');
    end
    allkeys = fieldnames(recon.param);
    if (size(recon.param.(allkeys{1}), 1) < labelmax) % if label-based, copy
        if (~isfield(recon, 'seg'))
            error('label-based param found in recon mesh, but recon.seg not defined');
        end
        cfg.param = recon.param; % need to be followed by rbupdateprop to propagate to cfg.prop
    else % node or element based param
        for i = 1:length(allkeys)
            cfg.param.(allkeys{i}) = meshinterp(recon.param.(allkeys{i}), recon.mapid, recon.mapweight, recon.elem, cfg.param.(allkeys{i}));
        end
        % need to be followed by rbupdateprop to propagate to cfg.prop
    end
elseif (isfield(recon, 'prop')) % map recon.prop to cfg.prop if param does not exist
    if (~isa(recon.prop, 'containers.Map')) % single wavelength
        if (size(recon.prop, 1) < labelmax) % if label-based, copy
            if (isfield(recon, 'node') && ~isfield(recon, 'seg'))
                error('label-based prop found in recon mesh, but recon.seg not defined');
            end
            cfg.prop = recon.prop;
        elseif (isfield(recon, 'mapid')) % if node/elem based, interpolate
            cfg.prop = meshinterp(recon.prop, recon.mapid, recon.mapweight, recon.elem, cfg.prop);
        end
    else % multiple wavelengths
        allkeys = recon.prop.keys;
        if (size(recon.prop(allkeys{1}), 1) < labelmax) % label based
            if (~isfield(recon, 'seg'))
                error('label-based param found in recon mesh, but recon.seg not defined');
            end
            cfg.prop = recon.prop;
        else % node/elem based
            for i = 1:length(allkeys)
                cfg.prop(allkeys{i}) = meshinterp(recon.prop(allkeys{i}), recon.mapid, recon.mapweight, recon.elem, cfg.prop(allkeys{i}));
            end
        end
    end
else
    warning('nothing to sync, recon.param and recon.prop properties are not found');
end
