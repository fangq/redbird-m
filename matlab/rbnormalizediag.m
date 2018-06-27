function [A, di]=rbnormalizediag(A0)

% normalize/equalize a square matrix and make the diagnonal all 1s

Adiag=diag(A0);
di=1./sqrt(Adiag);
A=(di*di').*A0;