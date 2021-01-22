function varargout=rbrun(cfg,recon,detphi0,varargin)
%
% detphi=rbrun(cfg)
% newrecon=rbrun(cfg,recon,detphi0)
% [newrecon,resid,newcfg,...]=rbrun(cfg,recon,detphi0,sd)
%   or
% [newrecon,resid,newcfg,...]=rbrun(cfg,recon,detphi0,sd,'param1',value1,'param2',value2,...)
%
% Master script to run streamlined forward, inversion with various
% reconstruction modes.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the forward data structure defining the forward mesh in a
%          reconstruction, if cfg is the only input, rbrun is the same as
%          rbrunforward
%     recon: if recon is a vector, it must be the measurement data for all
%               src/detector pairs
%          if recon is a struct, it can have the following subfields
%               recon.node, recon.elem: the reconstruction mesh
%               recon.lambda: Tikhonov regularization parameter, if not
%                     present, use 0.1
%               recon.prop or recon.param: initial guess of mua/mus or
%                     parameters
%     detphi0: measurement data to be fitted, can be either a matrix or a
%              containers.Map (multi-wavelength); if detphi0 is a struct
%              with subfields like node/elem, it is treated as a forward
%              simulation data structure like cfg, and is used to run the
%              forward solution to obtain the simulated "measurement"
%
% output:
%     newrecon,resid,newcfg...: if only cfg is used as input, the output
%              parameters are the same as in rbrunforward; if recon is
%              given, the output sequence is the same as that of rbrunrecon
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

mode='image';
opt=struct;
if(nargin==4 || ~ischar(varargin{1}) || bitand(length(varargin),1)==1)
    sd=varargin{1};
    if(length(varargin)>=3)
        opt=varargin2struct(varargin{2:end});
    end
else
    opt=varargin2struct(varargin{:});
end

mode=jsonopt('mode',mode,opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run forward for the heterogeneous domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin==1)
    [varargout{1:nargout}]=rbrunforward(cfg);
    return;
elseif(nargin==2)
    if(~isfield(recon,'detphi0'))
        error('you must give detphi0 as input');
    else
        detphi0=recon.detphi0;
    end
end
if(isstruct(detphi0) && isfield(detphi0,'srcpos'))
    detphi0=rbrunforward(detphi0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Reset the domain to a homogeneous medium for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch mode
    case {'bulk','seg'}
        if(strcmp(mode,'bulk'))
            recon.seg=ones(size(recon.node,1),1);
            maxseg=1;
        else
            maxseg=max(recon.seg);
        end
        if(isfield(recon,'param') && isstruct(recon.param) && isfield(recon,'bulk'))
            for ch=fieldnames(recon.param)
               recon.param.(ch{1})=recon.bulk.(ch{1})*ones(maxseg,1);
            end
            % cfg.prop will be updated inside rbrunrecon;
        end
        if((~isfield(recon,'param') && isfield(recon,'prop') && isempty(recon.prop)) || ~isfield(recon,'prop'))
            nref=1.37;
            if(isfield(recon.bulk,'n'))
                nref=recon.bulk.n;
            end
            if(isfield(recon.bulk,'musp'))
                recon.prop=repmat([recon.bulk.mua,recon.bulk.musp,0,nref],maxseg,1);
            else
                recon.prop=repmat([recon.bulk.mua,1/(3*recon.bulk.dcoeff),0,nref],maxseg,1);
            end
            recon.prop=[0 0 1 1; recon.prop];
        end
end

if(isfield(recon,'seg') && ~isempty(recon.seg))
    if(numel(recon.seg)==1 || length(unique(recon.seg))==1)
        recon.seg=ones(size(recon.node,1),1);
        fprintf('rbrun: run bulk estimation ...\n');
    elseif(isvector(recon.seg))
        fprintf('rbrun: run segmented tissue property estimation ...\n');
    elseif(min(size(recon.seg))>=2)
        if(strcmp(mode,'fuzzy'))
            fprintf('rbrun: run fuzzy-segmentation based property estimation ...\n');
        elseif(strcmp(mode,'prior'))
            fprintf('rbrun: run composition-prior guided reconstruction ...\n');
        elseif(strcmp(mode,'roi'))
            fprintf('rbrun: run roi-based composition-prior guided reconstruction ...\n');
        end
    end
end

if(~exist('sd','var'))
    sd=rbsdmap(cfg);
end

maxiter=jsonopt('maxiter',10,opt);

[varargout{1:nargout}]=rbrunrecon(maxiter,cfg,recon,detphi0,sd,opt);
