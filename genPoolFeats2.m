% Given a set of element detectors and an image, generate
% a feature vector suitable for use in an SVM.  Specifically,
% run the detectors on the image, and find the max detection
% in each cell of a 2-level spatial pyramid (2x2 and 1x1).  Hence
% the final representation has length numel(detr.id)*5.  All features
% for one detector are stored consecutively.
function [feature] = genPoolFeats(detr, im)
global ds;
    pyramid = constructFeaturePyramid(im, ds.conf.params);
    [features, levels, indexes,gradsums] = unentanglePyramid(pyramid, ...
      ds.conf.params.patchCanonicalSize/ds.conf.params.sBins-2);
  invalid=(gradsums<9);
  size(features)
  features(invalid,:)=[];
  levels(invalid)=[];
  indexes(invalid,:)=[];
  gradsums(invalid)=[];
  disp(['threw out ' num2str(sum(invalid)) ' patches']);
  patsz=ds.conf.params.patchCanonicalSize;
  fsz=(patsz-2*ds.conf.params.sBins)/ds.conf.params.sBins;
  pos=pyridx2pos(indexes,reshape(levels,[],1),fsz,pyramid);
  posy=(pos.y1 + pos.y2)/2+.000001;
  posx=(pos.x1 + pos.x2)/2+.000001;
  idx=1;
  feature=zeros(5,size(detr.w,1));
  % find the maximum match in each of the 4 grid cells separately
  % (note assigntoclosest finds the best score in a set of feature vectors
  % for each detector).
  for(i=[-1 1])
    for(j=[-1 1]);
      [~,dist]=assigntoclosest(detr.w,features(i*(posy-size(im,1)/2) > 0 & j*(posx-size(im,2)/2) > 0,:),1);
      dist=dist(:)-detr.b;
      feature(idx,:)=dist(:)';
      idx=idx+1;
    end
  end
  % compute the 1x1 layer of the spatial pyramid as the max
  % of the 2x2 layer.
  feature(end,:)=max(feature(1:end-1,:),[],1);
  feature=feature(:)';
end
