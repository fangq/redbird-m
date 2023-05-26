function com = wfCOM(wf,cfg)

node = cfg.node;

wf = wf';

x = sum(wf.*repmat(node(:,1),1,size(wf,2)),1);
y = sum(wf.*repmat(node(:,2),1,size(wf,2)),1);
z = sum(wf.*repmat(node(:,3),1,size(wf,2)),1);

com = [x(:) y(:) z(:)];