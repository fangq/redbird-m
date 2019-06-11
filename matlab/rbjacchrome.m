function Jchrome=rbjacchrome(Jmua, extin)

% convert the Jacobian for mua to chromorphores
% extin: rows: wavelengths; columns: chromorphores

extin=extin';
extin=extin(:)';
if(~ismatrix(Jmua))
    Jmua=rbmatflat(Jmua); %J=[J(wavelen1); J(wavelen2) ;...; J(wavelenN)]
end
Jchrome=kron(extin, Jmua);

% J=[J(chrom1),J(chrom2),...,J(chromM)]