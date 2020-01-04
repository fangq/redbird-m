function srcbc=rbsrc2bc(cfg,isdet)
%
% srcbc=rbsrc2bc(cfg)
%
% Converting wide-field source forms into a boundary condition by defining
% in-ward flux on the mesh surface
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the simulation data structure, with srctype, srcpos, srcdir, 
%          srcparam1, srcparam2, srcpattern fields
%     isdet: default is 0; if set to 1, rbsrc2bc process widefield
%          detectors, the relevant fields are dettype, detpos, detdir,
%          detparam1, detparam2, detpattern
%
% output:
%     srcbc: an array of Ns x Nt, where Nt is size(cfg.face,1) and Ns is
%            the number of sources (if isdet=1, the detector counts)
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details 
%
% -- this function is part of Redbird-m toolbox
%

srcbc=[];

if(nargin<2)
        isdet=0;
end

if(~isdet)
        if(~isfield(cfg,'srctype') || strcmp(cfg.srctype,'pencil') || strcmp(cfg.srctype,'isotropic'))
                return;
        end
        srctype=cfg.srctype;
        srcpos=cfg.srcpos;
        srcdir=cfg.srcdir;
        srcparam1=cfg.srcparam1;
        srcparam2=cfg.srcparam2;
        if(strcmp(srctype,'pattern'))
                srcpattern=cfg.srcpatten;
        end
else
        if(~isfield(cfg,'dettype') || strcmp(cfg.dettype,'pencil') || strcmp(cfg.dettype,'isotropic'))
                return;
        end
        srctype=cfg.dettype;
        srcpos=cfg.detpos;
        srcdir=cfg.detdir;
        srcparam1=cfg.detparam1;
        srcparam2=cfg.detparam2;
        if(strcmp(srctype,'pattern'))
                srcpattern=cfg.detpatten;
        end
end

% already converted
if(size(srcpos,2)==size(cfg.face,1))
        srcbc=srcpos;
        return;
end

switch srctype
    case {'planar','pattern','fourier'}
            ps=[srcpos; srcpos+srcparam1(1:3); ...
                srcpos+srcparam1(1:3)+srcparam2(1:3); srcpos+srcparam2(1:3); srcpos];
            c0=meshcentroid(cfg.node,cfg.face);
            newnode=rotatevec3d([c0; ps],srcdir(1:3));
            srcpoly=newnode(end-4:end, 1:2);
            [isin,ison]=inpolygon(newnode(1:end-5,1),newnode(1:end-5,2),srcpoly(:,1),srcpoly(:,2));
            isin=isin | ison;
            idx=find(isin);
            if(~isempty(idx)) % the below test only works for convex shapes
                    AB=cfg.node(cfg.face(idx,2),1:3)-cfg.node(cfg.face(idx,1),1:3);
                    AC=cfg.node(cfg.face(idx,3),1:3)-cfg.node(cfg.face(idx,1),1:3);
                    N=cross(AB',AC')';
                    dir=sum(N.*repmat(srcdir(:)',size(N,1),1),2);
                    if(all(dir>=0))
                            error('please reorient the surface triangles');
                    end
                    srcbc=zeros(1,size(cfg.face,1));
                    srcbc(idx(dir<0))=1;
                    pbc=newnode(idx(dir<0),1:2);
                    
                    dp=pbc-repmat(srcpoly(1,:),size(pbc,1),1);
                    dx=srcpoly(2,:)-srcpoly(1,:);
                    dy=srcpoly(4,:)-srcpoly(1,:);
                    nx=dx/norm(dx);
                    ny=dy/norm(dy);
                    
                    bary=[sum(dp.*repmat(nx/norm(dx),size(dp,1),1),2), sum(dp.*repmat(ny/norm(dy),size(dp,1),1),2)];
                    bary(bary<0)=0;
                    bary(bary>=1)=1-1e-6;
                    
                    if(exist('srcpattern','var'))
                            if(ismatrix(srcpattern))
                                    srcpattern=permute(srcpattern,[3 1 2]);
                            end
                            pdim=size(srcpattern);
                            patsize=pdim(1);
                            srcbc=repmat(srcbc,patsize,1);
                            for i=1:patsize
                                    srcbc(idx(dir<0),i)=srcpattern(sub2ind(pdim, i*ones(size(bary,1),1), floor(bary(:,1)*pdim(1)), floor(bary(:,2)*pdim(2))));
                            end
                    elseif(strcmp(srctype,'fourier'))
                            kx=floor(srcparam1(4));
                            ky=floor(srcparam2(4));
                            phi0=(srcparam1(4)-kx)*2*pi;
                            M=1-(srcparam2(4)-ky);
                            srcbc=repmat(srcbc,kx*ky,1);
                            for i=1:kx
                                    for j=1:ky
                                        srcbc((i-1)*ky+j,idx(dir<0))=0.5*(1+M*cos((i*bary(:,1)+j*bary(:,2))*2*pi+phi0))';
                                    end
                            end
                    end
            else
                    error('source direction does not intersect with the domain');
            end
    otherwise
            error('this source type is not supported');
end
