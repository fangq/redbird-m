function [Anew, allkeys]=rbmatflat(Amat)

allkeys=[];
if(ismatrix(Amat))
    Anew=Amat;
    return;
end

if(isa(Amat,'containers.Map'))
    allkeys=Amat.keys;
    for key=allkeys
        id=key{1};
        Anew=[Anew; Amat(id)];
    end
elseif(isstruct(Amat))
    allkeys=fieldnames(Amat);
    for key=allkeys
        id=key{1};
        Anew=[Anew; Amat.(id)];
    end
end