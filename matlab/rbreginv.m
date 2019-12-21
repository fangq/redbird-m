function res=rbreginv(Amat, varargin)
%
% res=rbreginv(Amat, rhs, lambda, Lqr)
%
% Solve a linear equation and automatically decide the least-square solution
% method (for underdetermined or overdetermined cases) based on the shape of
% the input data
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

if(size(Amat,1)>=size(Amat,2)) % overdetermined case
    res=rbreginvover(Amat, varargin{:});
else                         % underdetermined case
    res=rbreginvunder(Amat, varargin{:});
end