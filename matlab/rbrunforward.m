function [detval, phi, Amat, rhs, sflag]=rbrunforward(cfg,varargin)
%
% [detval, phi]=rbrunforward(cfg)
%    or
% [detval, phi, Amat, rhs]=rbrunforward(cfg,'param1',value1,...)
%
% Perform forward simulations at all sources and all wavelengths based on the input structure
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the redbird data structure
%
% output:
%     detval: the values at the detector locations
%     phi: the full volumetric forward solution computed at all wavelengths
%     Amat: the left-hand-side matrices (a containers.Map object) at specified wavelengths 
%     rhs: the right-hand-side vectors for all sources (independent of wavelengths)
%     param/value pairs: (optional) additional parameters
%          'solverflag': a cell array to be used as the optional parameters
%               for rbfemsolve (starting from parameter 'method'), for
%               example  rbrunforward(...,'solverflag',{'pcg',1e-10,200})
%               calls rbfemsolve(A,rhs,'pcg',1e-10,200) to solve forward
%               solutions
% 
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(~isfield(cfg,'deldotdel'))
    cfg.deldotdel=rbdeldotdel(cfg);
end

wavelengths={''};
if(isa(cfg.prop,'containers.Map'))
   wavelengths=cfg.prop.keys;
end

Amat=containers.Map();
phi=containers.Map();
detval=containers.Map();
opt=varargin2struct(varargin{:});
solverflag=jsonopt('solverflag',{},opt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rhs,loc,bary]=rbfemrhs(cfg);

for waveid=wavelengths
	wv=waveid{1};

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%   Build LHS
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Amat(wv)=rbfemlhs(cfg,cfg.deldotdel,wv); % use native matlab code, 1 sec for 50k nodes

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%   Solve for solutions at all nodes: Amat*res=rhs
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%solverflag={'pcg',1e-12,200}; % if iterative pcg method is used
	[phi(wv),sflag]=rbfemsolve(Amat(wv),rhs,solverflag{:});

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%   Extract detector readings from the solutions
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	detval(wv)=rbfemgetdet(phi(wv), cfg, loc, bary); % or detval=rbfemgetdet(phi(wv), cfg, rhs); 
end

% if only a single wavelength is required, return regular arrays instead of a map
if(length(wavelengths)==1)
    Amat=Amat(wavelengths{1});
    phi=phi(wavelengths{1});
    detval=detval(wavelengths{1});
end