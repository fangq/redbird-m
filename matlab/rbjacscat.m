function Jscat=rbjacscat(Jd, dcoeff, scatpow, wv)
%
% Jscat=rbjacchrome(Jmua, dcoeff, wv)
%
% Building the Jacobian matrices for scattering-amplitude and
% scattering-power from Jacobian of the diffusion coeff
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Jd: the Jacobian for diffusion coefficient, as a containers.Map
%     dcoeff: the current diffusion coeff at each node
%     scatpow: the current scattering power at each node
%     wv: wavelength list, if not given, Jd must be a containers.Map
%
% output:
%     Jscat: the Jacobian in a struct as Jscat.{scatamp,scatpow}
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(nargin<4)
    wv=keys(Jd);
end

% Jscat=[J(scatamp),J(scatpow)]
Jscat=struct('scatamp',[],'scatpow',[]);
% Jscat.scatamp=zeros(size(cell2mat(Jd.values'),1),size(cell2mat(Jd.values'),2));
% Jscat.scatpow=Jscat.scatamp;

for i=1:length(wv)
%     Jscat.scatamp(((i-1)*size(Jd(wv{i}),1)+1):(i)*size(Jd(wv{i}),1),:) =...
%         rbjacscatamp(Jd(wv{i}), dcoeff(wv{i}), str2double(wv{i}), scatpow);
%     Jscat.scatpow(((i-1)*size(Jd(wv{i}),1)+1):(i)*size(Jd(wv{i}),1),:) =...
%         rbjacscatpow(Jd(wv{i}), dcoeff(wv{i}), str2double(wv{i}));
    Jscat.scatamp = [Jscat.scatamp; rbjacscatamp(Jd(wv{i}), dcoeff(wv{i}), str2double(wv{i}), scatpow)];
    Jscat.scatpow = [Jscat.scatpow; rbjacscatpow(Jd(wv{i}), dcoeff(wv{i}), str2double(wv{i}))];
end
