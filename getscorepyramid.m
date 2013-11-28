function [res,assignedidx,dist,clustid]=getscorepyramid(w,features,indexes,thresh,invalid,conf)
  try
  %inprod=bsxfun(@plus,single(features)*single(w)',single(thresh(:)'));
  inprod=single(features)*single(w)';
  inprod(invalid,:)=-Inf;
  [y,x]=find(bsxfun(@gt,inprod,thresh(:)'));
  clustid=reshape(x,[],1);
  assignedidx=y;
  dist=reshape(inprod(sub2ind(size(inprod),y(:),x(:))),[],1);
  [inprod,indexes]=distributeby(inprod,indexes(:,1:2),indexes(:,3));
  for(i=1:numel(indexes))
    sz{i}=max(indexes{i},[],1);
    indexes{i}=sub2ind(sz{i},indexes{i}(:,1),indexes{i}(:,2));
  end
  for(i=size(inprod,1):-1:1)
    for(j=size(inprod{i},2):-1:1)
      res{j,1}{i,1}(indexes{i})=inprod{i}(:,j);
      res{j,1}{i,1}=reshape(res{j}{i},sz{i});
    end
  end
  catch ex,dsprinterr;end
end

