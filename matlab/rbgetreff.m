function Reff = rbgetreff(n_in,n_out)

% given refractive index of the diffuse medium, calculate the effective
% refractive index, defined as in Haskell.

% original file name calcExtBnd
% author: David Boas <dboas at bu.edu>
% this file was modified from the PMI toolbox
% License: BSD 3-clause license.

if(nargin==1)
    n_out=1;
end

oc = asin(1/n_in);
ostep = pi / 2000;

o = 0:ostep:oc;

cosop = (1-n_in^2 * sin(o).^2).^0.5;
coso = cos(o);
r_fres = 0.5 * ( (n_in*cosop-n_out*coso)./(n_in*cosop+n_out*coso) ).^2;
r_fres = r_fres + 0.5 * ( (n_in*coso-n_out*cosop)./(n_in*coso+n_out*cosop) ).^2;

r_fres(ceil(oc/ostep):1000) = 1;

o = 0:ostep:ostep*(length(r_fres)-1);
coso = cos(o);

r_phi_int = 2 * sin(o) .* coso .* r_fres;
r_phi = sum(r_phi_int) / 1000 * pi/2;

r_j_int = 3 * sin(o) .* coso.^2 .* r_fres;
r_j = sum(r_j_int) / 1000 * pi/2;

Reff = (r_phi + r_j) / (2 - r_phi + r_j);
