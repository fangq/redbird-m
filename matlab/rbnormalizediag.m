function [A, Adiag]=rbnormalizediag(A0)

% normalize/equalize a square matrix and make the diagnonal all 1s

Adiag=diag(A0);
Adiag=1./sqrt(Adiag);
A=(Adiag*Adiag').*A0;