function res=rbreginvover(Amat, rhs, lambda, Lqr, varargin)

% solve an overdetermined Gauss-Newton normal equation
%  delta_mu=inv(J'J + lambda*L'L)*J'*(y-phi)

rhs=Amat'*rhs(:);

Hess=Amat'*Amat; % Gauss-Hessian matrix, approximation to Hessian (2nd order)

Hess(1:1+size(Hess,1):end)=Hess(1:1+size(Hess,1):end)+lambda;

if(nargin<4 || isempty(Lqr))
    Hess(1:1+size(Hess,1):end)=Hess(1:1+size(Hess,1):end)+lambda;
else
    Hess=Hess+lambda*(Lqr'*Lqr);
end

res=rbfemsolve(Hess, rhs, varargin{:});
