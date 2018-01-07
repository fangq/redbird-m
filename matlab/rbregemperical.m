function lambda=rbregemperical(Hess, residual, alpha)

% an emperical approach to determine the regularization coeff lambda

ggav=mean(diag(Hess));
lambda=alpha*ggav*residual*residual;