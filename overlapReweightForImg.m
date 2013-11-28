% compute the alpha values for a single image.  dets is just a list 
% of detections in standard form [x1 y1 x2 y2 score detrid imageid flip]
function reweights=overlapReweightForImg(dets,detridclassgroup)
global ds;
[~,idx]=ismember(dets(:,6),detridclassgroup(:,1));
imgs=dsload('ds.imgs{ds.conf.currimset}');
toovlweight1=detridclassgroup(idx,2)==imgs.label(dets(:,7));
if(~dsbool(ds.conf.params,'ovlweight'))
  reweights=ones(size(dets(:,1)));
else
  toovlweightf=find(toovlweight1);
  % break the detections up into flipped and un-flipped.
  if(size(dets,2)>=8)
    toovlweightall{1}=toovlweightf(dets(toovlweightf,8)==0);
    toovlweightall{2}=toovlweightf(dets(toovlweightf,8)==1);
  else
    toovlweightall={toovlweightf};
  end
  reweights=ones(size(idx));

  for(i=1:numel(toovlweightall))
    toovlweight=toovlweightall{i};
    if(numel(toovlweight)>0)
      conf.groups=detridclassgroup(idx,3);
      conf.groups=conf.groups(toovlweight);
      % we set the initial weights of each detection to the detection score;
      % these weights in effect determine how strongly each detection "competes"
      % for pixels.
      % Detections with a negative score get a small epsilon value,
      % because we're computing the re-weigting based on the inverses of
      % these scores.  Hence negative values make no sense.
      % The output of overlapReweight is scaled by the input weight (overlapReqeight RE-weights
      % it doesn't know that we made up the input weights).
      reweights(toovlweight)=overlapReweight(dets(toovlweight,1:4),max(0,dets(toovlweight,5))+.000001,ds.conf.params.sBins,ds.conf.params.patchCanonicalSize,conf)./(max(0,dets(toovlweight,5))+.000001);
      if(any(reweights>1e10)),error('fail');end
    end
  end
end
% finally, scale the weights by the number of detections for each
% detector.
udetr=unique(idx);
idx(~toovlweight1)=0;
for(i=1:numel(udetr))
  reweights(idx==udetr(i))=reweights(idx==udetr(i))./sum(idx==udetr(i));
end
