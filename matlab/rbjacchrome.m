function Jchrome=rbjacchrome(Jmua, extin)

% convert the Jacobian for mua to chromorphores
% extin: rows: wavelengths; columns: chromorphores

extin=extin';
extin=extin(:)';
Jchrome=kron(extin, Jmua);

% J=[J(wavelen1,chrom1),J(wavelen1,chrom2),...,J(wavelen2,chrom1)...]