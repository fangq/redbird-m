function Amat=rbdefemlhs(cfg)

nn=size(cfg.node,1);
ne=size(cfg.elem,1);
nf=size(cfg.face,1);
ntot=length(cfg.cols);

Amat=sparse(nn,nn);

mua=cfg.prop(cfg.elemprop+1,1);
mus=cfg.prop(cfg.elemprop+1,2);
g=cfg.prop(cfg.elemprop+1,3);
musp=(1-g).*mus;

Dcoeff=1./(3.0*(mua+musp));

count=1;
for i=1:4
    for j=i:4
        ra=Dcoeff.*cfg.deldotdel(:,count).*cfg.evol*(1/120);
        if(j==i)
            Amat=Amat+ra;
        else
            Aoffd=Aoffd+ra;
        end
        count=count+1;
    end
end
    