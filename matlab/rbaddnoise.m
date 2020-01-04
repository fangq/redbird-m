function newdata=rbaddnoise(data, snrshot, snrthermal,randseed)
%
% newdata=rbaddnoise(data, snrshot, snrthermal,randseed)
%
% Adding simulated shot-noise and thermal noise to the simulated data
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     data: noise-less data
%     snrshot: the desired SNR of the shot-noise (amplitude-proportional)
%     snrthermal: the desired SNR of the thermal-noise (noise floor)
%     randseed (optional): specify a seed for reproducible noise
%
% output:
%     newdata: noise contaminated data
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

if(nargin<4)
    randseed=123456789;
end

if(nargin>=3)
    rand('state',randseed);
    randn('state',randseed);
end

if(nargin==1)
    newdata=data;
    warning('no noise added');
    return;
end

if(nargin<3)
    snrthermal=inf;
    if(nargin<2)
        snrshot=inf;
    end
end

datanorm=abs(data);
max_amp=max(datanorm(:));

sigma_shot=10^(-real(snrshot)/20);
sigma_thermal=max_amp*10^(-real(snrthermal)/20);

if(isreal(data))
    newdata=data + sqrt(abs(data)).*randn(size(data))*sigma_shot + randn(size(data))*sigma_thermal;
else
    sigma_shot_phase=10^(-image(snrshot)/20);
    sigma_thermal_phase=10^(-image(snrthermal)/20);
    ampshotnoise=randn(size(data))*sigma_shot;
    phaseshotnoise=randn(size(data))*sigma_shot_phase*2*pi;
    ampthermalnoise=randn(size(data))*sigma_thermal;
    phasethermalnoise=randn(size(data))*sigma_thermal_phase*2*pi;
    shotnoise=sqrt(abs(data)).*complex(ampshotnoise.*cos(phaseshotnoise), ampshotnoise.*sin(phaseshotnoise));
    thermalnoise=complex(ampthermalnoise.*cos(phasethermalnoise), ampthermalnoise.*sin(phasethermalnoise));
    newdata=data + shotnoise + thermalnoise;
end
    