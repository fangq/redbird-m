function [cfg, recon, resid, updates, Jmua, detphi, phi]=rbrunrecon(maxiter,cfg,recon,detphi0,sd,varargin)
%
% [cfg, recon, resid]=rbrunrecon(maxiter,cfg,recon,detphi0,sd)
%   or
% [cfg, recon, resid, updates, Jmua, detphi, phi]=rbrunrecon(maxiter,cfg,recon,detphi0,sd,'param1',value1,'param2',value2,...)
%
% Perform a single iteration of a Gauss-Newton reconstruction
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     maxiter: number of iterations
%     cfg: simulation settings stored as a redbird data structure
%     recon: reconstruction data structure, recon may have
%         node: reconstruction mesh node list
%         elem: reconstruction mesh elem list
%         param: wavelength-independent parameter on the recon mesh
%         prop: wavelength-dependent optical properties on the recon mesh
%         lambda: Tikhonov regularization parameter
%         mapid: the list of element indices of the reconstruction mesh where each forward mesh
%           node is enclosed
%         mapweight: the barycentric coordinates of the forward mesh nodes inside the
%           reconstruction mesh elements
%     detphi0: measurement data vector or matrix
%     sd (optional): source detector mapping table, if not provided, call
%         rbsdmap(cfg) to compute
%     param/value: acceptable optional parameters include
%         'reform': 'complex' (default),'logphase',or 'real'
%         'lambda': Tikhonov regularization parameter (0.05), overwrite recon.lambda
%         'report': 1 (default) to print residual and timing; 0 silent mode
%         'tol': convergence tolerance, if relative residual is less than
%                this value, stop, default is 0, which runs maxiter
%                iterations
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
    return;
end

resid=zeros(1,maxiter);
updates=repmat(struct,1,maxiter);

opt=varargin2struct(varargin{:});

lambda=0.05;
if(isfield(recon,'lambda'))
    lambda=recon.lambda;
end
lambda=jsonopt('lambda',lambda,opt);

doreport=jsonopt('report',1,opt);
convergetol=jsonopt('tol',0,opt);
reform=jsonopt('reform','complex',opt);

if(nargin<5)
    sd=rbsdmap(cfg);
end

% start iterative Gauss-Newton based reconstruction

for iter=1:maxiter
    tic

    % update forward mesh prop/param using recon mesh prop/param if given
    % for rbsyncprop to work, one must provide initial values of cfg.prop
    % (or cfg.param) if recon.prop (or recon.param) is specified
    if(isfield(recon,'node') && isfield(recon,'elem'))
        [cfg,recon]=rbsyncprop(cfg,recon);
    end

    % run forward on forward mesh, cfg.param, if present, is propagated to 
    % cfg.prop when calling rbfemlhs inside rbrunforward
    [detphi, phi]=rbrunforward(cfg);
    
    % build Jacobians on forward mesh
    if(isfield(cfg,'omega') && cfg.omega>0) % if RF data
        % currently, only support node-based values; rbjac supports
        % multiple wavelengths, in such case, it returns a containers.Map
        if((isfield(cfg,'seg') && length(cfg.seg)==size(cfg.elem,1)) || size(cfg.prop,1)==size(cfg.elem,1))
            % element based properties
            [Jmua_n, Jmua, Jd_n, Jd]=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol);
            clear Jmua_n Jd_n % do not use ~ because older Octave does not support
        else
            % node based properties
            [Jmua, Jmua_e, Jd]=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol);
            clear Jmua_e
        end
    else % CW only
        if((isfield(cfg,'seg') && length(cfg.seg)==size(cfg.elem,1)) || size(cfg.prop,1)==size(cfg.elem,1))
            [Jmua_n, Jmua]=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol);
            clear Jmua_n;
        else
            Jmua=rbjac(sd, phi, cfg.deldotdel, cfg.elem, cfg.evol);
        end
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
    if(isfield(recon,'elem') && isfield(recon,'node') && isfield(recon,'mapid') && isfield(recon,'mapweight')) % dual-mesh reconstruction
        Jmua=structfun(@(x) transpose(meshremap(x.',recon.mapid, recon.mapweight,recon.elem,size(recon.node,1))), Jmua,'UniformOutput',false); 
    end

    if(isfield(recon,'seg')) % reconstruction of segmented domains
        Jmua=structfun(@(x) rbmasksum(x,recon.seg(:)'), Jmua,'UniformOutput',false);
    end

    % blocks contains unknown names and Jacob size, should be Nsd*Nw rows
    % and Nn columns, Nsd is src/det pairs, Nw is number of wavelengths,
    % and Nn is the recon mesh, if present, or forward mesh node size.
    blocks=structfun(@size, Jmua, 'UniformOutput', false);
    
    % flatten Jmua into a horizontally/column contatenated matrix
    if(strcmp(reform,'complex')==0)
        [Jflat,misfit]=rbmatreform(rbmatflat(Jmua), detphi0(:), detphi(:), reform);
    else
        Jflat=rbmatflat(Jmua);
        misfit=detphi0(:)-detphi(:);
    end

    % store the residual
    resid(iter)=sum(abs(misfit));
    updates(iter).detphi=detphi;
   
    % solver the inversion (J*delta_x=delta_y) using regularized
    % Gauss-Newton normal equation
    dmu_recon=rbreginv(Jflat, misfit, lambda);  % solve the update on the recon mesh
   
    % obtain linear index of each output species
    len=cumsum([1; structfun(@(x) x(2), blocks)]);
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
                            meshinterp(recon.prop(:,propidx),recon.mapid, recon.mapweight,recon.elem,cfg.prop(:,propidx)); % interpolate the update to the forward mesh
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
                error('unknown type %s is not supported',output{i});
        end
    end
    if(doreport)
        fprintf(1,'iter [%4d]: residual=%e, relres=%e lambda=%e (time=%f s)\n',iter, resid(iter), resid(iter)/resid(1), lambda, toc);
    end
    if(iter>1 && abs((resid(iter)-resid(iter-1))/resid(1))<convergetol)
        resid=resid(1:iter);
        updates=updates(1:iter);
        break;
    end
end

if(isfield(recon,'node') && isfield(recon,'elem'))
    [cfg,recon]=rbsyncprop(cfg,recon);
    cfg.prop=rbupdateprop(cfg);
end