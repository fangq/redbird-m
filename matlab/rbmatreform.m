function [newA, newrhs, nblock]=rbmatreform(Amat, ymeas, ymodel, output)

% convert an inverse problem to the log-amplitude and unwrapped phase form
%  

if(nargin<4)
    output='complex';
end

wavelengths={''};

if(isa(Amat,'containers.Map'))
    wavelengths=Amat.keys;
else
    Amat=containers.Map({''},{Amat});
    ymeas=containers.Map({''},{ymeas});
    ymodel=containers.Map({''},{ymodel});
end

newA=containers.Map();
newrhs=containers.Map();

nblock=1;

if(strcmp(output,'complex'))
    newA=Amat;
    for wv=wavelengths
        newrhs(wv)=ymeas(wv)-ymodel(wv);
    end
else
    for wv=wavelengths
        rhs=ymeas(wv)-ymodel(wv);
        if(strcmp(output,'real'))
            newA(wv)=real(Amat(wv));
            newrhs(wv)=real(rhs);

            if(~isreal(rhs) && ~isreal(Amat(wv)))
                newA(wv)=[real(Amat(wv)) -imag(Amat(wv)); 
                      imag(Amat(wv)) real(Amat(wv))];
                newrhs(wv)=[real(rhs); imag(rhs)];
                nblock=2;
            end
        elseif(strcmp(output,'logphase'))
            temp=repmat(conj(ymodel(wv))./abs(ymodel(wv).*ymodel(wv)),1,size(Amat(wv),2)).*Amat(wv);
            if(isreal(ymodel(wv)))
                newA(wv)=real(temp);
                newrhs(wv)=log(abs(ymeas(wv)))-log(abs(ymodel(wv)));
            else
                newA(wv)=[real(temp) ; imag(temp)];
                newrhs(wv)=[log(abs(ymeas(wv)))-log(abs(ymodel(wv))); 
                            angle(ymeas(wv))   -angle(ymodel(wv))];
                nblock=2;
            end
        end
    end
end


% if only a single wavelength is required, return regular arrays instead of a map
if(length(wavelengths)==1)
    newA=newA(wavelengths{1});
    if(nargout>1)
        newrhs=newrhs(wavelengths{1});
    end
end