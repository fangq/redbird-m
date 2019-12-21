function [newA, newrhs, nblock]=rbmatreform(Amat, ymeas, ymodel, output)
%
% [newA, newrhs, nblock]=rbmatreform(Amat, ymeas, ymodel, output)
%
% Reformat the matrix equation A*x=(ymeas-ymodel) to choose between the log-amplitude/phase
% form or the complex form
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Amat: the LHS of the matrix equation
%     ymeas: the vector that stores the measured data for the model to fit 
%     ymodel: the model predicted measurements at all source detector pairs, 
%           with a length matching that of ymeas
%
% output:
%     newA: the reformed LHS matrix 
%     newrhs: the reformed RHS matrix
%     nblock: nblock=length(newrhs)/size(Amat,1)
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(nargin<4)
    output='complex';
end

rhs=ymeas-ymodel;

nblock=1;

if(strcmp(output,'complex'))
    newA=Amat;
    newrhs=rhs;
    return;
end

if(strcmp(output,'real'))
    newA=real(Amat);
    newrhs=real(rhs);

    if(~isreal(rhs) && ~isreal(Amat))
        newA=[real(Amat) -imag(Amat); 
              imag(Amat) real(Amat)];
	    newrhs=[newrhs; imag(rhs)];
        nblock=2;
    end
    return;
end

if(strcmp(output,'logphase'))
    temp=repmat(conj(ymodel)./abs(ymodel.*ymodel),1,size(Amat,2)).*Amat;
    if(isreal(ymodel))
        newA=real(temp);
        newrhs=log(abs(ymeas))-log(abs(ymodel));
    else
        newA=[real(temp) ; imag(temp)];
	    newrhs=[log(abs(ymeas))-log(abs(ymodel)); 
                angle(ymeas)   -angle(ymodel)];
        nblock=2;
    end
    return;
end
