function res=rbreginvover(Amat, rhs, lambda, Lqr, varargin)
%
% res=rbreginvover(Amat, rhs, lambda, Lqr, varargin)
%
% Solve an overdetermined Gauss-Newton normal equation
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Amat: the left-hand-side matrices (a containers.Map object) at specified wavelengths 
%     rhs: the right-hand-side vectors for all sources (independent of wavelengths)
%     lambda: the Tikhonov regularization parameter
%     Lqr: the QR decomposition of the L-matrix used as the regularization matrix
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
%  delta_mu=inv(J'J + lambda*L'L)*J'*(y-phi)

rhs=Amat'*rhs(:);

Hess=Amat'*Amat; % Gauss-Hessian matrix, approximation to Hessian (2nd order)

[Hess,Gdiag]=rbnormalizediag(Hess);

if(nargin<4 || isempty(Lqr))
    Hess(1:1+size(Hess,1):end)=Hess(1:1+size(Hess,1):end)+lambda;
else
    Hess=Hess+lambda*(Lqr'*Lqr);
end

res=Gdiag(:).*rbfemsolve(Hess, Gdiag(:).*rhs, varargin{:});
