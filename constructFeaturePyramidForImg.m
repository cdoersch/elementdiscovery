function pyramid = constructFeaturePyramidForImg(im, params, levels)
% levels: What level of pyramid to compute the features for.
%
% Author: saurabh.me@gmail.com (Saurabh Singh).
sBins = params.sBins;
'constructfeatpyr'
size(im)
I = im;
if(dsfield(params,'imageCanonicalSize'))
  canonicalSize = params.imageCanonicalSize;
  [IS, canoScale] = convertToCanonicalSize(I, canonicalSize);
else
  canonicalSize=min(size(I(:,:,1)));
  canoScale=1;
  IS=I;
end
size(IS)
[rows, cols, chans] = size(IS);
if chans < 3
  I = repmat(I, [1 1 3]);
  fprintf('WARNING: Image has < 3 channels, replicating channels\n');
end

numLevels  = getNumPyramidLevels(rows, cols, params.scaleIntervals, ...
  params.patchCanonicalSize)
scales = getLevelScales(numLevels, params.scaleIntervals);
if nargin < 3 || isempty(levels)
  levels = 1 : numLevels;
end
if(dsbool(params,'useColor'))
  im2=RGB2Lab(im).*.0025;
  if(dsbool(params,'extraColor'))
    im2=im2.*25;
  end
end
if(dsbool(params,'useColorHists'))
  im2=RGB2Lab(im);
end
if(dsbool(params,'patchOnly'))
  im2=rgb2gray(im);
end
pyramidLevs = cell(1, numLevels);
histbin=-100:20:100;
histbin(end)=histbin(end)+1;
gradientLevs={};
for i = 1 : length(levels)
  lev = levels(i);
  I1 = imresize(I, canoScale / scales(lev),'bilinear');
  [nrows, ncols, unused_dims] = size(I1);
  rowRem = rem(nrows, sBins);
  colRem = rem(ncols, sBins);
  if rowRem > 0 || colRem > 0
    I1 = I1(1:nrows-rowRem, 1:ncols-colRem, :);
  end
  if(~dsbool(params,'patchOnly'))
    if(dsbool(params,'overallnorm'))
      feat=features_nonorm(I1,sBins)*.02;
    else
      feat = features(I1, sBins);
    end
    [rows,cols,~]=size(feat);
  else
    rows=round(size(im2,1)/sBins);
    cols=round(size(im2,2)/sBins);
  end
  if(dsbool(params,'patchOnly'))
    feat=imresize(im2,[rows cols],'bilinear');
  else
    feat=feat(:,:,1:31);
    if(dsbool(params,'maxmin'))
      feat=cat(3,feat(:,:,1:18),max(feat(:,:,1:9),feat(:,:,10:18)),...
                                min(feat(:,:,1:9),feat(:,:,10:18)),feat(:,:,28:end));
    end
    if(dsbool(params,'medrat'))
       f=features_medrat(I1,sBins)+.000001;
       feat=cat(3,feat,f(2:end-1,2:end-1,:));
    end
    if(dsbool(params,'useColor'))
      feat=cat(3,feat,imresize(im2(:,:,2),[rows cols],'bilinear'),imresize(im2(:,:,3),[rows cols],'bilinear'));
    elseif(dsbool(params,'useColorHists'))
      for(k=2:3)
        tohist=im2col(imresize(im2(:,:,k),[rows*sBins, cols*sBins]),[sBins,sBins],'distinct');
        fhist{k-1}=reshape(permute(histc(tohist,histbin,1),[2,1]),rows,cols,[]);
        fhist{k-1}=fhist{k-1}(:,:,1:end-1).*.0006;
      end
      feat=cat(3,feat,fhist{1},fhist{2});
    end
  end
  [GX, GY] = gradient(I1);
  GI = mean((GX*255).^2, 3) + mean((GY*255).^2, 3);
  GI=imresize(GI,[rows,cols],'bilinear');
  pyramidLevs{lev} = feat;%(:, :, 1:31);
  gradientLevs{lev} = GI;
end
canoSize.nrows = size(im,1);
canoSize.ncols = size(im,2);
pyramid = struct('features', {pyramidLevs}, 'scales', scales, ...
  'canonicalScale', canoScale, 'sbins', sBins, 'canonicalSize', canoSize, 'gradimg', {gradientLevs});
end

function numLev = getNumPyramidLevels(rows, cols, intervals, basePatSize)
lev1 = floor(intervals * log2(rows / basePatSize(1)));
lev2 = floor(intervals * log2(cols / basePatSize(2)));
numLev = min(lev1, lev2) + 1;
end
