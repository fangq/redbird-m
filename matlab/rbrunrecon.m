function [cfg, recon, resid, updates, Jmua, detphi, phi]=rbrunrecon(maxiter,cfg,recon,detphi0,f2rid,f2rweight,reform)
%
% [cfg, recon, resid, Jmua]=rbrunrecon(maxiter,cfg,recon,detphi0,f2rid,f2rweight, reform)
%
% Perform a single iteration of a Gauss-Newton reconstruction
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     maxiter: number of iterations
%     cfg: simulation settings stored as a redbird data structure
%     sd: source detector mapping table
%     recon: reconstruction data structure, recon may have
%         node: reconstruction mesh node list
%         elem: reconstruction mesh elem list
%         param: wavelength-independent parameter on the recon mesh
%         prop: wavelength-dependent optical properties on the recon mesh
%         lambda: Tikhonov regularization parameter
%     detphi0: measurement data vector or matrix
%     f2rid: the list of element indices of the reconstruction mesh where each forward mesh
%           node is enclosed
%     f2rweight: the barycentric coordinates of the forward mesh nodes inside the
%           reconstruction mesh elements
%
% output:
%     cfg: the updated cfg structure, containing forward mesh and
%          reconstructed values in cfg.prop or cfg.param
%     recon: the updated recon structure, containing recon mesh and
%          reconstructed values in recon.prop or recon.param
%     resid: the residual betweet the model and the measurement data for
%          each iteration
%     Jmua: Jacobian in a struct form, each element is the Jacobian of an
%          unknown block
%     
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(maxiter==0 && nargin<3)
    % return detphi as cfg and phi as recon
    [cfg,recon]=rbrunforward(cfg);
end

resid=zeros(1,maxiter);
updates=repmat(struct,1,maxiter);

% start iterative Gauss-Newton based reconstruction

for iter=1:maxiter

    % update forward mesh prop/param using recon mesh prop/param if given
    if(isfield(recon,'node') && isfield(recon,'elem'))
        if(isfield(recon,'param')) % map recon.param to cfg.param
            cfg.param=structfun(@(x) meshinterp(x,f2rid, f2rweight,recon.elem), recon.param, 'UniformOutput', false);
        elseif(isfield(recon,'prop')) % map recon.prop to cfg.prop if param does not exist
            if(size(recon.prop,1)<min(size(recon.node,1),size(recon.elem,1))) % if label-based, copy
                cfg.prop=recon.prop;
            else % if node/elem based, interpolate
                cfg.prop=structfun(@(x) meshinterp(x,f2rid, f2rweight,recon.elem), recon.prop, 'UniformOutput', false);
            end
        end
    end

    % run forward on forward mesh, cfg.param, if present, is propagated to 
    % cfg.prop when calling rbfemlhs inside rbrunforward
    [detphi, phi]=rbrunforward(cfg);
    
    % build Jacobians on forward mesh
    if(isfield(cfg,'omega') && cfg.omega>0) % if RF data
        % currently, only support node-based values; rbjac supports
        % multiple wavelengths, in such case, it returns a containers.Map
        [Jmua, ~, Jd]=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol);
    else % CW only
        Jmua=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol);
    end
    % Jmua/Jd are either containers.Map(wavelength) or single matrix
    
    % build Jacobians for chromophores in the form of a struct
    % TODO: need to handle Jmua is a map but cfg.param is not defined
    if(isa(Jmua,'containers.Map') && isfield(cfg,'param') && isa(cfg.param,'struct'))
        if(exist('Jd','var'))
            [Jmua,detphi0,detphi]=rbmultispectral(Jmua, detphi0, detphi, cfg.param, Jd, cfg.prop);
        else
            [Jmua,detphi0,detphi]=rbmultispectral(Jmua, detphi0, detphi, cfg.param);
        end
    else  % recon mua/d per wavelengths
        if(exist('Jd','var'))
            Jmua=struct('mua',Jmua,'dcoeff',Jd);
        else
            Jmua=struct('mua',Jmua);
        end
    end

    % here, Jmua is a struct of unknown species; wavelengths
    % are vertically/row-wise concatenate; detphi and detphi0 are the
    % concatenated model and measurement RHS vectors

    % mapping jacobians from forward mesh to reconstruction mesh
    if(isfield(recon,'elem') && isfield(recon,'node') && nargin>6) % dual-mesh reconstruction
        Jmua=structfun(@(x) transpose(meshremap(x.',f2rid, f2rweight,recon.elem,size(recon.node,1))), Jmua); 
    end

    % blocks contains unknown names and Jacob size, should be Nsd*Nw rows
    % and Nn columns, Nsd is src/det pairs, Nw is number of wavelengths,
    % and Nn is the recon mesh, if present, or forward mesh node size.
    blocks=structfun(@size, Jmua, 'UniformOutput', false);
    
    % flatten Jmua into a horizontally/column contatenated matrix
    if(nargin>=8)
        [Jflat,misfit]=rbmatreform(rbmatflat(Jmua), detphi0(:), detphi(:), reform);
    else
        Jflat=rbmatflat(Jmua);
        misfit=detphi(:)-detphi0(:);
    end

    % store the residual
    resid(iter)=sum(abs(misfit));

    lambda=0.1;
    if(isfield(recon,'lambda'))
        lambda=recon.lambda;
    end
    
    % solver the inversion (J*delta_x=delta_y) using regularized
    % Gauss-Newton normal equation
    dmu_recon=rbreginv(Jflat, misfit, lambda);  % solve the update on the recon mesh
   
    % obtain linear index of each output species
    len=cumsum(structfun(@(x) x(2), blocks));
    output=fieldnames(blocks);
    for i=1:length(output)
        dx=dmu_recon(len(i):len(i+1)-1);
        updates(iter).(output{i})=dx;
        switch output{i}
            case {'mua','dcoeff'}
                propidx=strcmp(output{i},'dcoeff')+1;
                if(strcmp(output{i},'dcoeff')) % converting from d to musp
                    dx=1/dx*(1/3);
                end
                if(isfield(recon,'node')) % update recon mesh prop
                    if(length(dx)==size(recon.prop,1)-1) % label based prop
                        recon.prop(2:end,propidx)=recon.prop(2:end, propidx)+dx;
                        cfg.prop=recon.prop;
                    else % nodal or element based prop
                        recon.prop(:,propidx)=recon.prop(:, propidx)+dx;
                        cfg.prop(:,propidx)=...
                            meshinterp(recon.prop(:,propidx),f2rid, f2rweight,recon.elem); % interpolate the update to the forward mesh
                    end
                else % update forward mesh prop if single mesh
                    startid=(length(dx)==size(cfg.prop,1)-1)+1;
                    cfg.prop(startid:end,propidx)=cfg.prop(startid:end, propidx)+dx;
                end
            case {'hbo','hbr','water','lipid','scatamp','scatpow'}
                if(isfield(recon,'node')) % update recon mesh prop
                    recon.param.(output{i})=recon.param.(output{i})+dx;
                else
                    cfg.param.(output{i})=cfg.param.(output{i})+dx;
                end
            otherwise
                error(sprintf('unknown type %s is not supported',))
        end
    end
end
