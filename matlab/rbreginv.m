function res=rbreginv(Amat, varargin)
%
% res=rbreginv(Amat, rhs, lambda)
%   or
% res=rbreginv(Amat, rhs, lambda, Areg, blocks)
% res=rbreginv(Amat, rhs, lambda, Areg, blocks,...)
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
%     Areg (optional): the regularization matrix, use identity matrix if
%         not given, if Areg is given, it can be one of the two cases:
%         - if over-determined, it is the inversion of the unknown covariance
%             matrix Cx, where Cx=inv(L'L) 
%         - if under-determined, it is the inversion of the R matrix from
%             qr(L) where Cx=inv(L'L)
%         Areg is always square
%     blocks (optional): the dimensions of the 2D submatrices of Amat, needed if
%         size(Areg,1)==size(Areg,2) is not the same as size(Amat,2)
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