% given a set of element detectors, generate a feature vector 
% for an image that's suitable for classification in an SVM.
% Just performs genPoolFeats2 on the both flips of the image
% and averages the resulting feature vectors.
function feats=distGenPoolFeats(detrs,imidx);
  global ds;
  i=imidx;
  im=im2double(getimg(ds,i));
  feats=genPoolFeats2(detrs,im);
  im=im(:,end:-1:1,:);
  feats=(feats+genPoolFeats2(detrs,im))/2;
  feats=feats(:)';
end
