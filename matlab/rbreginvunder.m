function res=rbreginvunder(Amat, rhs, lambda, invLTL)

% solve an overdetermined Gauss-Newton normal equation
%  delta_mu=inv(L'L)*J'*inv(J*J' + lambda*I)*(y-phi)

rhs=rhs(:);

if(nargin>=4)
    Hess=Amat*invLTL*Amat'; % Gauss-Hessian matrix, approximation to Hessian (2nd order)
else
    Hess=Amat*Amat'; % Gauss-Hessian matrix, approximation to Hessian (2nd order)
end

Hess(1:1+size(Hess,1):end)=Hess(1:1+size(Hess,1):end)+lambda;

res=rbfemsolve(Hess, rhs);
res=Amat'*res;
if(nargin>=4)
    res=invLTL*res;
end
