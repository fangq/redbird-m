function [Anew, allkeys]=rbmatflat(Amat)
%
% [Anew, allkeys]=rbmatflat(Amat)
%
% Flatten a containers.Map object (in the case of multi-wavelength cases) into a rectangular 2D matrix
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Amat: the LHS of the matrix equation
%
% output:
%     Anew: the flattened LHS matrix 
%     allkeys: the wavelength labels in the Amat input if it is a containers.Map, or []
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

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