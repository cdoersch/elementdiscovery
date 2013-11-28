% Run detection on an image.  model is the standard detector format--
% i.e. it has a 'w' field where each row is a detector weight vector,
% a 'b' field with a bias that's subtracted from the score for each patch,
% and the 'id' field is an id for the detector.  the image is the image
% at index imid in .ds.imgs{ds.conf.currimset}.  conf can include:
%  - 'flipall': flip the images to get more detections. Default false
%  - 'thresh': detection threshold.  Default -Inf
%  - 'multperim': allow multiple detections per detector per image. Default false.
% If any of these fields aren't specified, they will be read from ds.conf.params
% before using the defaults.
%
% Output is a set of detections in the standard format: one per row,
% each row is [x1 y1 x2 y2 score detector_id image_id flip boxid].  Note
% that in this distribution, flip and boxid are not read.
function [dets,feats]=detectInIm(model,imid,conf)
  if(~exist('conf','var'))
    conf=struct();
  end
  if(~dsfield(conf,'thresh'))
    conf.thresh=-Inf;
  end
  conf2=conf;
  conf2.thresh=model.b+conf.thresh;
  boxid=[];
  [pos,dist,clustid,feats,flip,boxid]=bestInImbb(model.w,imid,conf2);
  dist=dist-model.b(clustid);
  if(isempty(dist))
    dets=zeros(0,9);
  else
    dets=[pos.x1,pos.y1,pos.x2,pos.y2,dist,model.id(clustid),repmat(imid,numel(dist),1),flip,boxid];
  end
end
