function varargout=rbplot(cfg, data, varargin)

if(size(data,2)==size(cfg.node,1) || size(data,2)==size(cfg.elem,1))
    [varargout{1:nargout}]=rbplotjacobian(cfg,data,varargin{:});
elseif(size(data,1)==size(cfg.node,1) || size(data,1)==size(cfg.elem,1))
    [varargout{1:nargout}]=rbplotforward(cfg,data,varargin{:});
else
    error('data row or column counts must match the node or elem number');
end