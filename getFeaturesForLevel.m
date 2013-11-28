function [features, indexes, gradsums] = getFeaturesForLevel(level, params,prSize, pcSize,pzSize,nExtra,gradimg)
% Author: saurabh.me@gmail.com (Saurabh Singh)
%edited by carl
[rows, cols, dims] = size(level);
%level=cat(3,level,imresize(im(:,:,2),[rows cols],'bilinear'),imresize(im(:,:,3),[rows cols],'bilinear'));
%dims=dims+2;
rLim = rows - prSize + 1;
cLim = cols - pcSize + 1;
%keyboard;
featDim = prSize * pcSize * pzSize + nExtra
if(dsbool(params,'appendbias'))
  features = ones(rLim * cLim, featDim+1);
else
  features = ones(rLim * cLim, featDim);
end
if(exist('gradimg','var'))
  gradsums=zeros(rLim * cLim,1);
  creategrad=true;
else
  creategrad=false;
end
%if(dsfield(params,'useColorHists'))
%  useColorHists=true;
%else
%  useColorHists=false;
%end
indexes = zeros(rLim * cLim, 2);
featInd = 0;
size(level)
for j = 1 : cLim
  for i = 1 : rLim
    feat = level(i:i+prSize-1, j:j+pcSize-1, :);
%    featlab = im2(i:i+prSize-1, j:j+pcSize-1, :)*.03;
    featInd = featInd + 1;
    %if(useColorHists)
    %  endhist=sum(sum(feat(:,:,32:end)));
    %  features(featInd,:)=[reshape(feat(:,:,1:31),1,[]) endhist(:)'];
    %else
      feat=reshape(feat, 1, []);
      %if(dsbool(params,'normalizefeats'))
      %  feat=(feat-mean(feat));
      %  feat=feat./max(.0000001,norm(feat));
      %end
      features(featInd, 1:size(feat,2)) = feat;

    %end
    if(creategrad)
      gradsums(featInd)=mean(mean(gradimg(i:i+prSize-1, j:j+pcSize-1)));
    end

%    features(featInd, :) = reshape(feat, 1, featDim);
    indexes(featInd, :) = [i, j];

  end
end
%if(dsbool(params,'normalizefeats'))
%  features=normrows(features,1,.000001);
%end
end
