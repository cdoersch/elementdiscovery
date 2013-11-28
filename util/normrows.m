function res=normrows(X,meansubtract,thresh)
  if(nargin>1&&meansubtract)
    X=bsxfun(@minus,X,mean(X,2));
  end
  if(nargin<3)
    thresh=0;
  end
  nrm=sqrt(sum(X.^2,2));
  badinds=nrm<=thresh;
  res=bsxfun(@rdivide,X,nrm);
  if(nargin>2)
    res(badinds,:)=0;
  end
end
