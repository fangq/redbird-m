function varargout=rbfemsolve(Amat, rhs, method, varargin)
if(nargin<3)
    method='mldivide';
end
if(nargout<1)
    error('output can not be empty');
end

% block solvers can handle multiple RHSs
switch method
    case 'blqmr'
        [varargout{:}]=blqmr(Amat,rhs,varargin{:});
        return;
    case 'mldivide'
        varargout{1}=Amat\rhs;
        return;
end


% non-block solvers have to solve one RHS at a time

sol=cell(1,nargout);
res=cell(1,nargout);

for i=1:size(rhs,2)
    switch method
        case 'qmr'
            [sol{:}]=qmr(Amat,rhs(:,i),varargin{:});
        case 'tfqmr'
            [sol{:}]=tfqmr(Amat,rhs(:,i),varargin{:});
        case 'cgs'
            [sol{:}]=cgs(Amat,rhs(:,i),varargin{:});
        case 'lsqr'
            [sol{:}]=lsqr(Amat,rhs(:,i),varargin{:});
        case 'gmres'
            [sol{:}]=gmres(Amat,rhs(:,i),varargin{:});
        case 'gmres'
            [sol{:}]=gmres(Amat,rhs(:,i),varargin{:});
        case 'pcg'
            [sol{:}]=pcg(Amat,rhs(:,i),varargin{:});
        case 'minres'
            [sol{:}]=minres(Amat,rhs(:,i),varargin{:});
        case 'symmlq'
            [sol{:}]=symmlq(Amat,rhs(:,i),varargin{:});
        case 'bicgstab'
            [sol{:}]=bicgstab(Amat,rhs(:,i),varargin{:});
        otherwise
            error(['method "',method,'" is not supported']);
    end
    for j=1:nargout
        res{j}=[res{j},sol{j}];
    end
end
varargout=res;