function [Jmua_n, Jmua_e, Jd_n, Jd_e]=rbjac(sd, phi, deldotdel, felem, evol)
%
% [Jmua_n, Jmua_e, Jd_n, Jd_e]=rbjac(sd, phi, deldotdel, felem, evol)
%
% Building the Jacobian matrices using adjoint method and native matlab code
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     sd: the source-detector mapping table
%     phi: the forward solutions at all sources
%     deldotdel: grad*phi dot product with grad phi, computed as part of the computation
%     felem: forward mesh element list
%     evol: forward mesh element volume
%
% output:
%     Jmua_n: the nodal Jacobian for absorption coeff. mua
%     Jmua_e: the element-wise Jacobian for absorption coeff. mua
%     Jd_n: (optional) the nodal Jacobian for diffusion coeff D
%     Jd_e: (optional) the element-wise Jacobian for diffusion coeff D
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(nargin<5 || isempty(sd) || isempty(phi) || isempty(deldotdel)|| isempty(evol))
    error('you must give at least 5 inputs and they must not be empty');
end

nelem=size(felem,1);

wavelengths={''};

if(isa(phi,'containers.Map'))
    wavelengths=phi.keys;
else
    phi=containers.Map({''},{phi});
end

Jmua_n=containers.Map();
Jmua_e=containers.Map();
if(nargout>2)
    Jd_n=containers.Map();
    Jd_e=containers.Map();
end

idx=[1 1 1 2 2 3 2 3 4 3 4 4
     2 3 4 3 4 4 1 1 1 2 2 3
     2 3 4 6 7 9 2 3 4 6 7 9];

for waveid=wavelengths
    wv=waveid{1};
    phiwv=phi(wv);
    if(isa(sd,'containers.Map'))
        sdwv=sd(wv);
    else
        sdwv=sd;
    end
    Jmua_node=zeros(size(sdwv,1),size(phiwv,1));
    Jmua_elem=zeros(size(sdwv,1),nelem);
    if(nargout>2)
        Jd_node=zeros(size(sdwv,1),size(phiwv,1));
        Jd_elem=zeros(size(sdwv,1),nelem);
    end
    for i=1:nelem
        phidotphi1=phiwv(felem(i,1:4),sdwv(:,1)).*phiwv(felem(i,1:4),sdwv(:,2));
        phidotphi2=phiwv(felem(i,idx(1,:)),sdwv(:,1)).*phiwv(felem(i,idx(2,:)),sdwv(:,2));
        Jmua_elem(:,i)=-(sum(phidotphi1,1)+sum(phidotphi2,1)*0.5)*(0.1*evol(i));
        if(nargout>2)
            Jcol=phidotphi1.'*deldotdel(i,[1 5 8 10]).';
            Jcol=Jcol+phidotphi2.'*deldotdel(i,idx(3,:)).';
            Jd_elem(:,i)=-Jcol;
        end
        for j=1:4
            Jmua_node(:,felem(i,j))=Jmua_node(:,felem(i,j))+Jmua_elem(:,i);
        end
        if(nargout>2)
            for j=1:4
                Jd_node(:,felem(i,j))=Jd_node(:,felem(i,j))-Jcol;
            end
        end
    end
    Jmua_node=Jmua_node*0.25;
    Jmua_n(wv)=Jmua_node;
    Jmua_e(wv)=Jmua_elem;
    if(nargout>2)
        Jd_node=Jd_node*0.25;
        Jd_n(wv)=Jd_node;
        Jd_e(wv)=Jd_elem;
    end
end

% if only a single wavelength is required, return regular arrays instead of a map
if(length(wavelengths)==1)
    Jmua_n=Jmua_n(wavelengths{1});
    Jmua_e=Jmua_e(wavelengths{1});
    if(nargout>2)
        Jd_n=Jd_n(wavelengths{1});
        Jd_e=Jd_e(wavelengths{1});
    end
end