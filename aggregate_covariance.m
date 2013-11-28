% Carl Doersch (cdoersch at cs dot cmu dot edu)
% compute features for a set of images, and then the sufficient
% statistics needed for computing the mean and covariance matrix.
ntotal=0;
for(i=1:numel(dsidx))
  params=ds.conf.params;
  I = im2double(getimg(ds,ds.myiminds(dsidx(i))));

  if(dsfield(params,'imageCanonicalSize'))
    [IS, scale] = convertToCanonicalSize(I, params.imageCanonicalSize);
  else
    IS=I;
  end
  pyramid = constructFeaturePyramidForImg(I, params);
  pcs=round(ds.conf.params.patchCanonicalSize/ds.conf.params.sBins)-2;
  [features, levels, indexes] = unentanglePyramid(pyramid, pcs,struct('normalizefeats',false));
  if(isempty(features))
    continue;
  end
  if(~exist('featsum','var'))
    featsum=sum(features,1);
  else
    featsum=featsum+sum(features,1);
  end
  ntotal=ntotal+size(features,1);
  if(~exist('dotsum','var'))
    dotsum=features'*features;
  else
    dotsum=dotsum+features'*features;
  end
end
if(ntotal==0)
  return
end
ds.n{dsidx(1)}=ntotal;
ds.featsum{dsidx(1)}=featsum;
ds.dotsum{dsidx(1)}=dotsum;
for(i=1:numel(dsidx))
  ds.imgflags{dsidx(i)}=1;
end
