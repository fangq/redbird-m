function Jmua=rbjacmua(sd, phi, nvol)

if(nargin<3 || isempty(sd) || isempty(phi) || isempty(nvol))
    error('you must give at least the first 3 inputs and they must not be empty');
end

if(size(phi,1)~=length(nvol))
    error('the row number of phi must be the same as the length of nvol');
end

Jmua=zeros(size(sd,1),size(phi,1));
for i=1:size(sd,1)
    Jmua(i,:)=phi(:,sd(i,1)).*phi(:,sd(i,2)).*nvol(:);
end