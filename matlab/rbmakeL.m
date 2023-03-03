function Lmat = rbmakeL (cfg,recon,prior,alpha,beta)

if (nargin < 4)
    alpha = 0.2;
end
if (nargin < 5)
    beta = 1.1;
end

tic
reconprior = [];
if (size(prior,1) == size(cfg.node,1))
%     [f2rid,f2rweight] = tsearchn(recon.node,recon.elem,cfg.node);
    reconprior = meshremap(prior,recon.mapid,recon.mapweight,recon.elem,size(recon.node,1));    
    reconprior = reconprior./sum(reconprior,2);
    reconprior(isnan(reconprior)) = 0;
    
    prior = reconprior;    
    clear reconprior
end


Lmat = zeros(size(prior,1),size(prior,1));

for ii = 1:size(prior,1)
    for jj = 1:size(prior,1)
        Lmat(ii,jj) = sum(abs(prior(ii,:) - prior(jj,:)));        % L1 Norm? Shouldn't this be L2?
%         Lmat(ii,jj) = sqrt(sum((prior(ii,:) - prior(jj,:)).^2)); 
        Lmat(jj,ii) = Lmat(ii,jj);
        
    end
%     fprintf(['looping ' num2str(ii) ' - ' num2str(toc) '\n']);       %   Debugging
end

Lmat(Lmat.^2 < ((alpha*size(prior,2)).^2)) = -alpha - Lmat(Lmat.^2 < ((alpha*size(prior,2)).^2))./size(prior,2);
Lmat(Lmat >= alpha*size(prior,2)) = 0;

Lmat(1:size(Lmat,2)+1:end) = 0;

% dmat = sqrt(abs(sum(Lmat,1)));
% dmat(dmat < 1e-5) = 1e-5;
% dmat = 1./dmat;
% Lmat = Lmat .* (1/beta) .* repmat(dmat,size(prior,1),1) .*repmat(dmat',1,size(prior,1));

dd = abs(sum(Lmat));
Lmat = Lmat./(beta.*sqrt(dd'*dd));

Lmat(1:size(Lmat,2)+1:end) = 1;


fprintf('Computing L matrix - %f seconds\n',toc);
end