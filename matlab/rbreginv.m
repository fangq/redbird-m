function res=rbreginv(Amat, varargin)

if(size(Amat,1)>=size(Amat,2)) % overdetermined case
    res=rbreginvover(Amat, varargin{:});
else                         % underdetermined case
    res=rbreginvunder(Amat, varargin{:});
end