function res=patch2feat(trpatches,conf)
  global ds;
  extrparams=ds.conf.params;
  extrparams.imageCanonicalSize=[max(min(size(trpatches{1}(:,:,1))),min(extrparams.patchCanonicalSize))];
  extrparams.basePatchSize=size(trpatches{1}(:,:,1));
  tmp=constructFeaturePyramidForImg(im2double(trpatches{1}),extrparams,1);
  %res(1,:)=tmp.features{1}(:)';
  pcs=round(ds.conf.params.patchCanonicalSize/ds.conf.params.sBins)-2;
  [res(1,:), levels, indexes,gradsums] = unentanglePyramid(tmp, ...
    pcs,conf);
  res(2:numel(trpatches),:)=0;%=zeros(numel(trpatches),numel(tmp.features{1}));
  for(i=2:numel(trpatches))
    tmp=constructFeaturePyramidForImg(im2double(trpatches{i}),extrparams,1);
    %res(i,:)=tmp.features{1}(:)';
    [res(i,:), levels, indexes,gradsums] = unentanglePyramid(tmp, ...
	pcs,conf);
  end
end
