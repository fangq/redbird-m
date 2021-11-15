function [newA, newrhs, nblock]=rbmatreform(Amat, ymeas, ymodel, form)
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
%     form: a string to indicate which form of output system, can be one of
%         'complex': no transformation, Amat and RHS keep original forms,
%             if inputs are real(complex), output are real(or complex)
%         'real': newA and newrhs are real matrices, x is assumed to be
%             real
%         'reim': newA and newrhs are real matrices, x is assumed to be
%             complex and is expanded to [real(x), imag(x)]' and newA is
%             expanded accordingly
%         'logphase': newA and newrhs are real matrices, x is assumed to be
%             complex and is expanded to [real(x) imag(x)]'; newrhs is
%             expanded as [log(ymeas-ymodel), angle(ymeans)-angle(ymodel)]'
%             and newA is transformed accordingly, please see Dr. Fang's
%             PhD Thesis, Section 3.2, Eq. 3.25
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
    form='complex';
end
nblock=1;

if (~isstruct(Amat))
    Amat = struct('J',Amat);
    ymeas = struct('detphi',ymeas);
    ymodel = struct('detphi',ymodel);
end

rfcw = length(Amat);
Aout = [];
rhsout = [];

for ii = rfcw:-1:1
    y0 = ymeas(ii).detphi;
    y1 = ymodel(ii).detphi;
    A = Amat(ii).J;
    
    rhs=y0-y1;

    if(strcmp(form,'complex'))
        newA=A;
        newrhs=rhs;
    end

    if(strcmp(form,'real') || strcmp(form,'reim'))
        newA=real(A);
        newrhs=real(rhs);

        if(~isreal(rhs) && ~isreal(A))
            if(strcmp(form,'reim'))
               newA=[real(A) -imag(A); 
                     imag(A) real(A)];
            else
               newA=[newA; imag(A)];
            end
            newrhs=[newrhs; imag(rhs)];
            nblock=1;
        end
    end

    if(strcmp(form,'logphase'))
        temp=repmat(conj(y1)./abs(y1.*y1),1,size(A,2)).*A;
        if(isreal(y1))
            newA=real(temp);
            newrhs=log(abs(y0))-log(abs(y1));
        else
            newA=[real(temp) ; imag(temp)];
            newrhs=[log(abs(y0))-log(abs(y1)); angle(y0)-angle(y1)];
            nblock=1;
        end
    end
    
    Aout = [Aout; newA];
    rhsout = [rhsout; newrhs];
end

newA = Aout;
newrhs = rhsout;
