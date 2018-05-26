function varargout=rbplotjacobian(cfg, jac, idx, varargin)

if(nargin<3)
    idx=1;
end

len=length(idx);

hd=zeros(len,1);

for i=1:len
    subplot(1,len,i);
    if(size(jac,2)==size(cfg.node,1))
        hd(i)=plotmesh([cfg.node, jac(idx(i),:)'], cfg.elem,varargin{:});
    elseif(size(jac,2)==size(cfg.elem,1))
        hd(i)=plotmesh(cfg.node, [cfg.elem(:,1:4), jac(idx(i),:)'],varargin{:});
    end
end

if(nargout>0)
    varargout{1}=hd;
end