function varargout=rbfemsolve(Amat, rhs, method, varargin)
%
% res=rbfemsolve(Amat, rhs, method, options ...)
%
% Solving a linear system defined by Amat*res=rhs using various methods
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the initial simulation data structure
%
% output:
%     Amat: the left-hand-size matrix, can be sparse
%     rhs:  the right-hand-size vector or matrix (multiple rhs)
%     method: (optional) a string specifying the solving method
%          mldivide: use the left-divide method: res=Amat\rhs
%          blqmr: use the block-QMR iterative method, other methods include
%          qmr, tfqmr, cgs, gmres, pcg, cgs, minres, symmlq, bicgstab
%     options: (optional) additional parameters are passed to the solver
%          function to specify solver options, please help qmr, cgs, etc to
%          see the additionally accepted parameters
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(nargin<3)
    method='mldivide';
end
if(nargout<1)
    error('output can not be empty');
end

% block solvers can handle multiple RHSs
switch method
    case 'blqmr'
        [varargout{1:nargout}]=blqmr(Amat,full(rhs),varargin{:});
        return;
    case 'mldivide'
        varargout{1}=full(Amat\rhs);
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
        res{j}=[res{j},full(sol{j})];
    end
end
varargout=res;