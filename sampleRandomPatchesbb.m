%v1 by Saurabh Singh
%Edit: Carl Doersch (cdoersch at cs dot cmu dot edu).
function [patches, patFeats, probabilities] = sampleRandomPatches(pos,samplelimit,conf)
if(~exist('conf','var'))
  conf=struct();
end
global ds;
params=ds.conf.params;
%data = pos;
%pos = pos.annotation;
rand('seed',1000*pos);
I2 = im2double(getimg(ds,pos));%imread([ds.conf.gbz{ds.conf.currimset}.cutoutdir ds.imgs(pos).fullpath]));
if(dsfield(conf,'detsforclass'))
  annot=getannot(pos);
  bbs=[annot.x1 annot.y1 annot.x2 annot.y2];
  classes=[annot.label];
  occl=annot.occluded;
  difficult=annot.difficult;
  boxid=annot.boxid;
  bbminsize=[bbs(:,4)-bbs(:,2)+1,bbs(:,3)-bbs(:,1)+1];
  valid=(~occl & c(ismember(classes,conf.detsforclass)) & ~difficult & all(bsxfun(@ge,bbminsize,ds.conf.params.patchCanonicalSize*1.5),2));
  bbs(~valid,:)=[];
  boxid(~valid)=[];
  bbs=[bbs,repmat([-Inf,-Inf,pos,0],size(bbs,1),1),boxid];
else
  bbs=[1,1,size(I2,2),size(I2,1),-Inf,-Inf,pos,0,1];
end
patchesall={};
patFeatsall={}
probabilitiesall={}
for(bbidx=1:size(bbs,1))
  I=I2(bbs(bbidx,2):bbs(bbidx,4),bbs(bbidx,1):bbs(bbidx,3),:);
  if(dsfield(params,'imageCanonicalSize'))
    [IS, scale] = convertToCanonicalSize(I, params.imageCanonicalSize);
  else
    IS=I;
    scale=1;
  end
  [rows, cols, unused] = size(IS);
  IG = getGradientImage(IS);
  pyramid = constructFeaturePyramidForImg(I, params);
  conf.imid=pos;
  pcs=round(ds.conf.params.patchCanonicalSize/ds.conf.params.sBins)-2;
  [features, levels, indexes] = unentanglePyramid(pyramid, ...
    pcs,conf);
    %params.patchCanonicalSize);
  selLevels = 1 : params.scaleIntervals/2 : length(pyramid.scales);
  levelScales = pyramid.scales(selLevels);
  numLevels = length(selLevels);

  patches = [];
  patFeats = [];
  probabilities = [];
  basenperlev=pyramid.features{selLevels(end)};
  basenperlev=(basenperlev(1)-pcs(1)+1)*(basenperlev(2)-pcs(2)+1);
  for i = 1 : numLevels
    levPatSize = floor(params.patchCanonicalSize .* levelScales(i));
    if(dsbool(params,'sampleBig'))
      numLevPat=floor(basenperlev/2);
    else
      numLevPat = floor((rows / levPatSize(1)) * ...
        (cols / levPatSize(2))*8);
    end
    %disp([num2str(levelScales(i)) '->' num2str(numLevPat)]);
    
    levelPatInds = find(levels == selLevels(i));
    if numLevPat <= 0
      continue;
    end
    
    IGS = IG;
    pDist = getProbDistribution(IGS, levPatSize);
    pDist1d = pDist(:);
    randNums = getRandForPdf(pDist1d, numLevPat);
    probs = pDist1d(randNums);
    [IY, IX] = ind2sub(size(IGS), randNums);
    IY = ceil(IY ./ (levelScales(i) * params.sBins));
    IX = ceil(IX ./ (levelScales(i) * params.sBins));
    
    [nrows, ncols, unused] = size(pyramid.features{selLevels(i)});
    IY = IY - floor(pcs(1) / 2);
    IX = IX - floor(pcs(2) / 2);
    xyToSel = IY>0 & IY<=nrows-pcs(1)+1 & IX>0 & IX<=ncols-pcs(2)+1;
    IY = IY(xyToSel);
    IX = IX(xyToSel);
    probs = probs(xyToSel);
    inds = sub2ind([nrows-pcs(1)+1 ncols-pcs(2)+1], IY, IX);
    [inds, m, unused] = unique(inds);
    probs = probs(m);
    selectedPatInds = levelPatInds(inds,:);
    fsz=(ds.conf.params.patchCanonicalSize-2*ds.conf.params.sBins)/ds.conf.params.sBins;
    metadata = pyridx2pos(indexes(selectedPatInds,:),levels(selectedPatInds),fsz,pyramid)
    %getMetadataForPositives(selectedPatInds, levels,...
    %  indexes, pcs(1), pcs(2), pyramid, pos);
    feats = features(selectedPatInds, :);
    if ~isempty(metadata)
      patInds = cleanUpOverlappingPatches(metadata, ...
      params.samplingOverlapThreshold);
      metadata.x1=metadata.x1+bbs(bbidx,1)-1;
      metadata.x2=metadata.x2+bbs(bbidx,1)-1;
      metadata.y1=metadata.y1+bbs(bbidx,2)-1;
      metadata.y2=metadata.y2+bbs(bbidx,2)-1;
      if(any(metadata.y2<metadata.y1))
        error('malformed patch');
      end
      n=size(metadata.x1);
      metadata.boxid=repmat(bbs(bbidx,9),n,1);
      metadata.flip=repmat(bbs(bbidx,8),n,1);
      metadata=effstridx(metadata,patInds);
      patches = [patches;[metadata.x1, metadata.y1, metadata.x2, metadata.y2, repmat([0 0],numel(metadata.x1),1), repmat(pos,numel(metadata.x1),1), metadata.flip, metadata.boxid]];
      patFeats = [patFeats; feats(patInds, :)];
      probabilities = [probabilities probs(patInds)'];
    end
    
  end
  patchesall{bbidx,1}=patches;
  patFeatsall{bbidx,1}=patFeats;
  probabilitiesall{bbidx}=probabilities;
end
patches=structcell2mat(patchesall);
patFeats=structcell2mat(patFeatsall);
probabilities=structcell2mat(probabilitiesall);
if(exist('samplelimit','var'))
  inds=randperm(size(patches,1));
  inds=inds(1:min(numel(inds),samplelimit));
  patches=patches(inds,:);
  patFeats=patFeats(inds,:);
  probabilities=probabilities(inds);
end
end

function patInds = cleanUpOverlappingPatches(patches, thresh)
  patmat=[patches.x1 patches.y1 patches.x2 patches.y2 rand(size(patches.y2))];
  patInds=myNms(patmat,thresh);
end

function [centers, vertExt] = getCategoryCenters(data, category)
objects = data.annotation.object;
objNames = {objects.name};
[ismem, unused] = ismember(objNames, {category});
primLoc = find(ismem);
centers = zeros(length(primLoc), 2);
vertExt = zeros(length(primLoc), 1);
for j = 1 : length(primLoc)
  vertExt(j) = getVerticalExtent(objects(primLoc(j)));
  [centers(j, 1), centers(j, 2)] = getCenter(objects(primLoc(j)), data);
end
end

function ext = getVerticalExtent(obj)
[x,y] = getLMpolygon(obj.polygon);
ext = max(y) - min(y) + 1;
end

function [cx cy] = getCenter(obj, data)
bb = getBoundingBox(obj, data.annotation);
cx = (bb(1) + bb(3)) / 2;
cy = (bb(2) + bb(4)) / 2;
end

function I1 = getGradientImage(I)
[GX, GY] = gradient(I);
I1 = sum(abs(GX), 3) + sum(abs(GY), 3);
I1 = I1.^2;
end

function dist = getProbDistribution(I, pSize)
h = fspecial('gaussian', pSize, min(pSize)/3);
I = imfilter(I, h);
dist = I ./ sum(sum(I));
end
