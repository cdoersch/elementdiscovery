%generally, targets are the detectors and toassign is patches from current image
function [assignedidx, dist, clustid]=findmatches(targets,toassign,thresh,conf);
  global ds;
  thresh=thresh(:);
  if(~exist('conf','var'))
    conf=struct();
  end
  if(dsfield(conf,'rho'))
    pertargrho=conf.rho(:);
  else
    pertargrho=0;
  end
  closest=zeros(size(toassign,1),1);
  outdist=zeros(size(toassign,1),1);
  assignedidx={};clustid={};dist={};
  for(i=1:800:size(toassign,1))
    inds=i:min(i+800-1,size(toassign,1));
    batch=toassign(inds,:);
    inprod=targets*(batch');
    inprod=bsxfun(@plus,inprod,pertargrho);
    [x,y]=find(bsxfun(@gt,inprod,thresh));
    clustid{end+1}=reshape(x,[],1);
    assignedidx{end+1}=reshape(inds(y),[],1);
    dist{end+1}=reshape(inprod(sub2ind(size(inprod),x(:),y(:))),[],1);
    %[outdist(inds),closest(inds)]=max(dist,[],1);
  end
  assignedidx=structcell2mat(assignedidx');
  clustid=structcell2mat(clustid');
  dist=structcell2mat(dist');
end
