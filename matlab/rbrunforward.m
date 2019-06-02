function [detval, phi, Amat, rhs]=rbrunforward(cfg)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Build RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rhs,loc,bary]=rbfemrhs(cfg);

for wv=wavelengths

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%   Build LHS
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Amat(wv)=rbfemlhs(cfg,cfg.deldotdel,wv); % use native matlab code, 1 sec for 50k nodes

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%   Solve for solutions at all nodes: Amat*res=rhs
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%phi(wv)=rbfemsolve(Amat(wv),rhs,'symmlq',1e-20,100);
	phi(wv)=rbfemsolve(Amat(wv),rhs);

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