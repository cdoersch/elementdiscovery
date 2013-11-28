function [features, levels, indexes, gradsums] = unentanglePyramid(pyramid, ...
  pcsz,conf)
  if(~exist('conf','var'))
    conf=struct();
  end
  global ds;
  conf=overrideConf(ds.conf.params,conf);
% Converts a pyramid of hog features for an image to a single matrix.
%
% Author: saurabh.me@gmail.com (Saurabh Singh).

  %prSize = round(patchCanonicalSize(1) / pyramid.sbins) - 2;
  %pcSize = round(patchCanonicalSize(2) / pyramid.sbins) - 2;
  %[prSize, pcSize, pzSize, nExtra]=getCanonicalPatchHOGSize(ds.conf.params);
  %prSize=pc
  selFeatures = cell(length(pyramid.features), 1);
  selGradsums = cell(length(pyramid.features), 1);
  selFeaturesInds = cell(length(pyramid.features), 1);

  selLevel = cell(length(pyramid.features), 1);
  totalProcessed = 0;
  if(numel(pyramid.features)==0)
    features=[];
    levels=[];
    indexes=[];
    gradsums=[];
    return
  end
  %im=RGB2Lab(im);
  if(numel(pcsz)<3)
    pcsz(3)=size(pyramid.features{1},3);
  end
  if(numel(pcsz)<4)
    pcsz(4)=0;
  end

  for i = 1 : length(pyramid.features)
    [feats, indexes, selGradsums{i,1}] = getFeaturesForLevel(pyramid.features{i}, conf, pcsz(1), ...
      pcsz(2),pcsz(3),pcsz(4),pyramid.gradimg{i});
    selFeatures{i} = feats;
    selFeaturesInds{i} = indexes;
    numFeats = size(feats, 1);
    selLevel{i} = ones(numFeats, 1) * i;
    totalProcessed = totalProcessed + numFeats;
  end
  gradsums=cell2mat(selGradsums);

  [features, levels, indexes] = appendAllTogether(totalProcessed, ...
    selFeatures, selLevel, selFeaturesInds);
  'featsize'
  size(feats)
  clear selFeatures;
  if(dsbool(conf,'normbeforewhit'))
    features=bsxfun(@rdivide,bsxfun(@minus,features,mean(features,2)),max(sqrt(var(features,1,2).*size(features,2)),.0000001));
  end
  if(dsbool(conf,'whitening'))
    try
      whit=dsload('.ds.whitenmat');
      whiten=1;
    catch
      whiten=0;
    end
    if(whiten)
      features=bsxfun(@minus,features,(dsload('.ds.datamean')))*(dsload('.ds.whitenmat')');
    end
  end
  if(dsbool(conf,'normalizefeats'))
    features=bsxfun(@rdivide,bsxfun(@minus,features,mean(features,2)),sqrt(var(features,1,2).*size(features,2)));
  end
  if(dsfield(conf,'contextfeats'))
    disp('contextfeats');
    dets=getImgDetections(conf.imid);
    if(~isempty(dets))
      dets=dets(ismember(dets(:,6),conf.contextfeats),:);
    end
    patsz=ds.conf.params.patchCanonicalSize;%allsz(resinds(k),:);
    fsz=(patsz-2*ds.conf.params.sBins)/ds.conf.params.sBins;
    %sz=size(ds.centers{c});
    %sz=sz(1:2);
    %idxpad=floor((sz-fsz)./2);
    imgs=getimgs();
    %pos=pyridx2pos(indexes,pyramid.canonicalScale,reshape(levels,[],1),...
    %       fsz(1),fsz(2),pyramid.sbins,[pyramid.canonicalSize.nrows pyramid.canonicalSize.ncols]);
    pos=pyridx2pos(indexes,reshape(levels,[],1),fsz,pyramid);
    pos=[pos.x1 pos.y1 pos.x2 pos.y2];
    features=[features contextfeats(pos,dets,conf.contextfeats)];
  end

end

function [newFeat, newLev, newInds] = appendAllTogether(totalProcessed, ...
  features, levels, indexes)
newFeat = zeros(totalProcessed, size(features{1}, 2));
newLev = zeros(totalProcessed, 1);
newInds = zeros(totalProcessed, 2);

featInd = 1;
for i = 1 : length(features)
  if isempty(features{i})
    continue;
  end
  startInd = featInd;
  endInd = startInd + size(features{i}, 1) - 1;
  newFeat(startInd:endInd, :) = features{i};
  features{i} = [];
  newLev(startInd:endInd) = levels{i};
  levels{i} = [];
  newInds(startInd:endInd, :) = indexes{i};
  indexes{i} = [];
  featInd = endInd + 1;
end
end
