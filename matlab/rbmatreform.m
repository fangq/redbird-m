function [newA, newrhs, nblock]=rbmatreform(Amat, ymeas, ymodel, output)

% convert an inverse problem to the log-amplitude and unwrapped phase form
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
