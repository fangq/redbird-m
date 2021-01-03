function res=rbreginvunder(Amat, rhs, lambda, invLTL)
%
% res=rbreginvunder(Amat, rhs, lambda, invLTL)
%
% Perform forward simulations at all sources and all wavelengths based on the input structure

% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Amat: the left-hand-side matrices (a containers.Map object) at specified wavelengths 
%     rhs: the right-hand-side vectors for all sources (independent of wavelengths)
%     lambda: the Tikhonov regularization parameter
%     invLTL: the inversion of the L-matrix used as regularization inv(L^TL)
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
%  delta_mu=inv(L'L)*J'*inv(J*J' + lambda*I)*(y-phi)

len=size(Amat,2);
idx=find(sum(Amat)~=0);
if(length(idx)<size(Amat,2))
    Amat=Amat(:,idx);
    %TODO: need to shrink invLTL as well
end

emptydata=find(sum(Amat')~=0);
if(length(emptydata)<size(Amat,1))
    Amat=Amat(emptydata,:);
    rhs=rhs(emptydata);
end

rhs=rhs(:);

if(nargin>=4)
    Hess=Amat*invLTL*Amat'; % Gauss-Hessian matrix, approximation to Hessian (2nd order)
else
    Hess=Amat*Amat'; % Gauss-Hessian matrix, approximation to Hessian (2nd order)
end

[Hess,Gdiag]=rbnormalizediag(Hess);

Hess(1:1+size(Hess,1):end)=Hess(1:1+size(Hess,1):end)+lambda;

res=Gdiag(:).*rbfemsolve(Hess, Gdiag(:).*rhs);

res=Amat'*res;
if(nargin>=4)
    res=invLTL*res;
end

if(length(idx)<len)
    res0=zeros(len,size(res,2));
    res0(idx,:)=res;
    res=res0;
end