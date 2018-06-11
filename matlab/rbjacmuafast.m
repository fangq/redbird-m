function Jmua=rbjacmuafast(sd, phi, nvol, elem)

if(nargin<3 || isempty(sd) || isempty(phi) || isempty(nvol))
    error('you must give at least the first 3 inputs and they must not be empty');
end

if(size(phi,1)==length(nvol))
    Jmua=zeros(size(sd,1),size(phi,1));
    for i=1:size(sd,1)
        Jmua(i,:)=phi(:,sd(i,1)).*phi(:,sd(i,2)).*nvol(:);
    end
elseif(nargin>3 && size(nvol,1)==size(elem,1))
    Jmua=zeros(size(elem,1),size(sd,1));
    for i=1:size(sd,1)
        for j=1:4
           Jmua(:,i)=Jmua(:,i)+phi(elem(:,j),sd(i,1)).*phi(elem(:,j),sd(i,2)).*nvol(:);
        end
    end
    Jmua=Jmua.'*0.25;
else
    error('the row number of phi must be the same as the length of nvol or elem is needed');
end

Jmua=-Jmua; % increasing mua, decreasing phi