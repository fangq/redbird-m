function [finalAmat, finalrhs, nblock]=rbcreateinv(Amat, ymeas, ymodel, params, output)

% convert an inverse problem to the log-amplitude and unwrapped phase form
%  

if(nargin<5)
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
    for waveid=wavelengths
        wv=waveid{1};
        newrhs(wv)=ymeas(wv)-ymodel(wv);
    end
else
    for waveid=wavelengths
        wv=waveid{1};
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
            temp=conj(ymodel(wv))./abs(ymodel(wv).*ymodel(wv));
            temp=repmat(temp(:),1,size(Amat(wv),2)).*Amat(wv);
            if(isreal(ymodel(wv)))
                newA(wv)=real(temp);
                temp=log(abs(ymeas(wv)))-log(abs(ymodel(wv)));
                newrhs(wv)=temp(:);
            else
                newA(wv)=[real(temp) ; imag(temp)];
                temp=log(abs(ymeas(wv)))-log(abs(ymodel(wv)));
                ptemp=angle(ymeas(wv))   -angle(ymodel(wv));
                newrhs(wv)=[temp(:); ptemp(:)];
                nblock=2;
            end
        end
    end
end

finalAmat=[];
finalrhs=[];

if(isa(params,'struct') && isfield())
    paramlist=fieldnames(params);
    extins=rbextinction(wavelengths,paramlist);
    for i=1:length(wavelengths)
        wv=waveid{i};
        wvAmat=[];
        for j=1:length(paramlist)
            wvAmat=[wvAmat, newA(wv)*extins(i,j)];
        end
        finalAmat=[finalAmat; wvAmat];
        finalrhs =[finalrhs; newrhs(wv)];
    end
else
    for i=1:length(wavelengths)
        wv=waveid{i};
        finalAmat=[finalAmat; newA(wv)];
        finalrhs =[finalrhs; newrhs(wv)];
    end
end