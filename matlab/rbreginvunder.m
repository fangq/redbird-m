function res=rbreginvunder(Amat, rhs, lambda, Lmat)

rhs=rhs(:);

Hess=Amat*Amat'; % Gauss-Hessian matrix, approximation to Hessian (2nd order)

Hess(1:1+size(Hess,1):end)=Hess(1:1+size(Hess,1):end)+lambda;

res=rbfemsolve(Hess, rhs);
res=LTL*Amat'*res;
