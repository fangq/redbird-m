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
                srcpattern=cfg.srcpattern;
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
                srcpattern=cfg.detpattern;
        end
end

% already converted
if(size(srcpos,2)==size(cfg.face,1))
        srcbc=srcpos;
        return;
end

z0=1/(cfg.prop(2,1)+cfg.prop(2,2)*(1-cfg.prop(2,3)));

switch srctype
    case {'planar','pattern','fourier'}
            ps=[srcpos; srcpos+srcparam1(1:3); ...
                srcpos+srcparam1(1:3)+srcparam2(1:3); srcpos+srcparam2(1:3); srcpos];
            
            pnode=cfg.node;
            pface=cfg.face;
            
            % if src is colimated (default), sink it by 1/mus'
            if(~isfield(cfg,'iscolimated') || cfg.iscolimated) 
                sinkplane=cfg.srcdir; % the plane where sinked planar source is located as [A,B,C,D] where A*x+B*y+C*z+D=0
                sinkplane(4)=-sum(cfg.srcdir.*(cfg.srcpos+cfg.srcdir*z0));
                [cutpos,cutvalue,facedata,elemid,nodeid]=qmeshcut(cfg.elem,cfg.node,zeros(size(cfg.node,1),1),sinkplane);
                pnode=cutpos;
                idx=find(facedata(:,3)~=facedata(:,4));
                pface=facedata(facedata(:,3)==facedata(:,4),1:3);
                pface=[pface;facedata(idx,[1 2 3]);facedata(idx,[1 3 4])];
            end
            c0=meshcentroid(pnode,pface);
            newnode=rotatevec3d([c0; ps],srcdir(1:3));
            srcpoly=newnode(end-4:end, 1:2);
            [isin,ison]=inpolygon(newnode(1:end-5,1),newnode(1:end-5,2),srcpoly(:,1),srcpoly(:,2));
            isin=isin | ison;
            idx=find(isin);
            if(~isempty(idx)) % the below test only works for convex shapes
                    AB=pnode(pface(idx,2),1:3)-pnode(pface(idx,1),1:3);
                    AC=pnode(pface(idx,3),1:3)-pnode(pface(idx,1),1:3);
                    N=cross(AB',AC')';
                    dir=sum(N.*repmat(srcdir(:)',size(N,1),1),2);
                    if(exist('sinkplane','var'))
                        dir(dir>0)=-dir(dir>0);
                    end
                    if(all(dir>=0))
                            error('please reorient the surface triangles');
                    end
                    srcbc=zeros(1,size(pface,1));
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
                            srcpattern=permute(srcpattern,[3 1 2]);
                            pdim=size(srcpattern);
                            patsize=pdim(1);
                            srcbc=repmat(srcbc,patsize,1);
                            for i=1:patsize
                                    srcbc(i,idx(dir<0))=srcpattern(sub2ind(pdim, i*ones(size(bary,1),1), floor(bary(:,1)*pdim(2))+1, floor(bary(:,2)*pdim(3))+1));
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

% at this point, srcbc stores the J- at each surface triangle (or sinked triangles)

Reff=cfg.reff;
maxbcnode=max(pface(:));
if(exist('nodeid','var'))
    nodeweight=nodeid(:,3);
    nodeid=nodeid(:,1:2);
    maxbcnode=max(max(nodeid(pface,:)));
    parea=elemvolume(pnode,pface);
else
    parea=cfg.area;
end

% 1/18 = 1/2*1/9, where 2 comes from the 1/2 in ls=(1+Reff)/(1-Reff)/2*D,
% and 1/9 = (1/6+1/12+1/12)/3, where A/6 is <phi_i,phi_j> when i=j, and
% A/12 is i!=j
Adiagbc=parea(:)*((1-Reff)/(18*(1+Reff)));
Adiagbc=repmat(Adiagbc,1,size(srcbc,1)).*(srcbc');
rhs=sparse(size(cfg.node,1),size(srcbc,1));

for i=1:size(srcbc,1)
    if(exist('nodeid','var'))
        allnodes=nodeid(pface,:);
        rhs(1:maxbcnode,i)=sparse(allnodes(:), 1, [Adiagbc(:,i).*nodeweight;Adiagbc(:,i).*(1-nodeweight)]);
    else
        rhs(1:maxbcnode,i)=sparse(cfg.face(:), 1, repmat(Adiagbc(:,i),1,3));
    end
    rhs(:,i)=rhs(:,i)/sum(rhs(:,i));
end
srcbc=rhs;
