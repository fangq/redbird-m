function res=rbreginvunder(Amat, rhs, lambda, invR, blocks, varargin)
%
% res=rbreginvunder(Amat, rhs, lambda, invR)
%
% Perform forward simulations at all sources and all wavelengths based on the input structure

% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Amat: the left-hand-side matrices (a containers.Map object) at specified wavelengths 
%     rhs: the right-hand-side vectors for all sources (independent of wavelengths)
%     lambda: the Tikhonov regularization parameter
%     invR: the inversion of the upper triangular matrix of QR-decomposed L
%
% output:
%     res: the least-square solution of the matrix equation
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

% solve an overdetermined Gauss-Newton normal equation
%  delta_mu=inv(L'L)*J'*inv(J*(inv(L'L))J' + lambda*I)*(y-phi)
%
% see Deng2015, to accelerate the computaiton, we do not take the inversion
% of L'L, instead, we perform QR decomposition of L, so that L=Q*R, and
% pass invR=inv(R) matrix here, this allows us to compute
%    Z=J*invR
%    J*(inv(L'L))J' = Z*Z' -> Hess
%    inv(L'L)*J'    = invR*Z'
%
% so we have
%    delta_mu=(invR*Z')*inv(Z*Z' + lambda*I)*(y-phi)

if(nargin>=4 && ~isempty(invR))
    nx=size(invR,1);
    if(nx==size(Amat,2))
        Amat=Amat*invR;  %Z=J*invR
    else
        if(isempty(blocks)) % assume the Hess matrix size is multiples of LTL
            for i=1:nx:size(Amat,1)
                Amat(:,i:i+nx-1)=Amat(:,i:i+nx-1)*invR;
            end
        else
            len=cumsum([1; structfun(@(x) x(2), blocks)]);
            for i=1:length(blocks)
                if(nx==len(i+1)-len(i)+1)  % if the block size match LTL
                    Amat(:,len(i):len(i+1)-1)=Amat(:,len(i):len(i+1)-1)*invR;
                end
            end
        end
    end
end

len=size(Amat,2);
idx=find(sum(Amat)~=0);
if(length(idx)<size(Amat,2))
    Amat=Amat(:,idx);
end

emptydata=find(sum(Amat')~=0);
if(length(emptydata)<size(Amat,1))
    Amat=Amat(emptydata,:);
    rhs=rhs(emptydata);
end

rhs=rhs(:);

Hess=Amat*Amat'; % Gauss-Hessian matrix, approximation to Hessian (2nd order)

[Hess,Gdiag]=rbnormalizediag(Hess);

Hess(1:1+size(Hess,1):end)=Hess(1:1+size(Hess,1):end)+lambda;

res=Gdiag(:).*rbfemsolve(Hess, Gdiag(:).*rhs, varargin{:});

if(nargin>=4 && ~isempty(invR))
    nx=size(invR,1);
    if(nx==size(Amat,2))
        res=invR*(Amat'*res);
    else
        if(isempty(blocks)) % assume the Hess matrix size is multiples of LTL
            for i=1:nx:size(Amat,1)
                res(i:i+nx-1)=invR*(Amat(:,i:i+nx-1)'*res(i:i+nx-1));
            end
        else
            len=cumsum([1; structfun(@(x) x(2), blocks)]);
            for i=1:length(blocks)
                if(nx==len(i+1)-len(i)+1)  % if the block size match LTL
                    res(len(i):len(i+1)-1)=invR*(Amat(:,len(i):len(i+1)-1)'*res(len(i):len(i+1)-1));
                end
            end
        end
    end
else
    res=Amat'*res;
end

if(length(idx)<len)
    res0=zeros(len,size(res,2));
    res0(idx,:)=res;
    res=res0;
end