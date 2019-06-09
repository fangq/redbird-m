function Jmua=rbjacmuafast(sd, phi, nvol, elem)

if(nargin<3 || isempty(sd) || isempty(phi) || isempty(nvol))
    error('you must give at least the first 3 inputs and they must not be empty');
end

wavelengths={''};

if(isa(phi,'containers.Map'))
    wavelengths=phi.keys;
else
    phi=containers.Map({''},{phi});
end

Jmua=containers.Map();

for waveid=wavelengths
    wv=waveid{1};
    phiwv=phi(wv);
    if(isa(sd,'containers.Map'))
        sdwv=sd(wv);
    else
        sdwv=sd;
    end
    if(size(phiwv,1)==length(nvol))
        Ja=zeros(size(sdwv,1),size(phiwv,1));
        for i=1:size(sdwv,1)
            Ja(i,:)=phiwv(:,sdwv(i,1)).*phiwv(:,sdwv(i,2)).*nvol(:);
        end
    elseif(nargin>3 && size(nvol,1)==size(elem,1))
        Ja=zeros(size(elem,1),size(sdwv,1));
        for i=1:size(sdwv,1)
            for j=1:4
               Ja(:,i)=Ja(:,i)+phiwv(elem(:,j),sdwv(i,1)).*phiwv(elem(:,j),sdwv(i,2)).*nvol(:);
            end
        end
        Ja=Ja.'*0.25;
    else
        error('the row number of phi must be the same as the length of nvol or elem is needed');
    end
    Jmua(wv)=-Ja; % increasing mua, decreasing phi
end

% if only a single wavelength is required, return regular arrays instead of a map
if(length(wavelengths)==1)
    Jmua=Jmua(wavelengths{1});
end
