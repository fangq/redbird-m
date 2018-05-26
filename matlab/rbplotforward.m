function varargout=rbplotforward(cfg, phi, idx, varargin)

if(nargin<3)
    idx=1;
end

len=length(idx);

hd=zeros(len,1);

for i=1:len
    subplot(1,len,i);
    hd(i)=plotmesh([cfg.node, phi(:,idx(i))], cfg.elem,varargin{:});
end

if(nargout>0)
    varargout{1}=hd;
end